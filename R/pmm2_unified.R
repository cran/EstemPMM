#' Unified PMM2 Estimator for Nonlinear Regression Models
#'
#' This file contains a unified PMM2 implementation suitable for arbitrary
#' nonlinear regression models (including SARIMAX) for which residuals
#' and their derivatives (Jacobian) can be computed.
#'
#' Two approaches are implemented:
#' 1. Iterative: Full iterative procedure (Nonlinear PMM2).
#' 2. One-step (Global): Single-step correction after a classical estimator (e.g., MLE).


#' Compute numerical Jacobian of the residual function
#'
#' Uses numDeriv::jacobian to compute the derivative matrix.
#'
#' @param fn_residuals Function(theta) returning the residual vector
#' @param theta Current parameter values
#' @param method Numerical differentiation method ("Richardson", "simple")
#' @return Jacobian matrix (n x p)
#' @keywords internal
compute_numerical_jacobian <- function(fn_residuals, theta, method = "Richardson") {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
        stop("Package 'numDeriv' is required for numerical Jacobian. Please install it.", call. = FALSE)
    }

    # Compute Jacobian: d(residuals)/d(theta)
    # For regression e = y - f(theta), so d(e)/d(theta) = -d(f)/d(theta)
    # numDeriv::jacobian computes d(fn)/d(theta)
    J_residuals <- numDeriv::jacobian(fn_residuals, theta, method = method)

    # For the PMM2 solver we need J = d(f)/d(theta) = -d(e)/d(theta)
    # Therefore we flip the sign
    J <- -J_residuals

    return(J)
}


#' Compute PMM2 weights and components
#'
#' @param residuals Residual vector
#' @return List with moments and PMM2 parameters (gamma3, gamma4, weights)
compute_pmm2_components <- function(residuals) {
    n <- length(residuals)

    # Central moments
    m1 <- mean(residuals)
    m2 <- mean((residuals - m1)^2)
    m3 <- mean((residuals - m1)^3)
    m4 <- mean((residuals - m1)^4)

    if (m2 < 1e-10) {
        return(NULL)
    } # Degenerate case

    # Standardized moments
    gamma3 <- m3 / m2^(1.5)
    gamma4 <- m4 / m2^2 - 3

    # PMM2 denominator check
    denom <- 2 + gamma4
    if (abs(denom) < 1e-6) denom <- 1e-6

    list(
        m2 = m2,
        gamma3 = gamma3,
        gamma4 = gamma4,
        denom = denom
    )
}

#' PMM2 step solver
#'
#' Solves the linearized system to find the parameter update.
#' Based on the Taylor expansion: e(theta) ~ e(theta_k) - J * delta
#'
#' @param residuals Current residuals
#' @param J Jacobian matrix (n x p), where J_ij = d(y_hat_i)/d(theta_j).
#'          We assume the standard regression definition: y = f(theta) + e,
#'          so e = y - f(theta) and d(e)/d(theta) = -d(f)/d(theta).
#'          We expect J = d(f)/d(theta) (gradient of the regression function).
#' @param pmm_stats Statistics from compute_pmm2_components
#'
#' @return Update vector delta
solve_pmm2_step <- function(residuals, J, pmm_stats) {
    # PMM2 maximizes (gamma3)^2 / (2 + gamma4).
    # The gradient of this function with respect to theta leads to the
    # estimating equations.
    #
    # Update formula (simplified for PMM2):
    # delta = (J'J)^(-1) J' * (residuals + correction)
    # where correction depends on the distribution asymmetry.
    #
    # Optimal estimating equation: sum( w(e_t) * grad(f_t) ) = 0
    # where w(e) is a polynomial weight function. For PMM2 (degree 2):
    # w(z) = z + a * (z^2 - 1)
    #
    # Coefficient 'a' is defined as:
    # a = - gamma3 / (2 + gamma4)
    #
    # Using the scoring approach where E[1 + 2*lambda*z] = 1:
    # delta = (J'J)^(-1) J' * (residuals + sigma * lambda * (z^2 - 1))

    sigma <- sqrt(pmm_stats$m2)
    z <- residuals / sigma

    # Coefficient for the quadratic term
    lambda <- -pmm_stats$gamma3 / pmm_stats$denom

    # Correction vector
    correction <- sigma * lambda * (z^2 - 1)

    # Effective "observations" for the LS step
    y_effective <- residuals + correction

    # LS step: delta = (J'J)^(-1) J' y_effective
    # Using qr.solve for numerical stability
    delta <- tryCatch(
        {
            qr.solve(J, y_effective)
        },
        error = function(e) {
            # Fallback if singular
            solve(t(J) %*% J + diag(1e-6, ncol(J))) %*% t(J) %*% y_effective
        }
    )

    as.vector(delta)
}


#' Universal PMM2 estimator (Iterative)
#'
#' @param theta_init Initial parameter values
#' @param fn_residuals Function(theta) returning the residual vector
#' @param fn_jacobian Function(theta) returning the Jacobian matrix (n x p).
#'                    J_ij = d(y_hat_i)/d(theta_j) = -d(epsilon_i)/d(theta_j).
#'                    If NULL, numerical Jacobian via numDeriv is used
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#'
#' @return List with results (theta, residuals, convergence, etc.)
#' @export
pmm2_nonlinear_iterative <- function(theta_init, fn_residuals, fn_jacobian = NULL,
                                     max_iter = 100, tol = 1e-6, verbose = FALSE) {
    theta <- theta_init
    p <- length(theta)

    if (verbose) cat("Starting Iterative PMM2...\n")

    for (iter in 1:max_iter) {
        # 1. Compute residuals and Jacobian at current point
        res <- fn_residuals(theta)

        if (is.null(fn_jacobian)) {
            J <- compute_numerical_jacobian(fn_residuals, theta)
        } else {
            J <- fn_jacobian(theta)
        }

        # Dimension check
        if (length(res) != nrow(J)) stop("Mismatch between residuals length and Jacobian rows")
        if (length(theta) != ncol(J)) stop("Mismatch between theta length and Jacobian cols")

        # 2. PMM2 statistics
        stats <- compute_pmm2_components(res)
        if (is.null(stats)) {
            warning("PMM2 stats computation failed (degenerate residuals).")
            break
        }

        # 3. Compute update step
        delta <- solve_pmm2_step(res, J, stats)

        # 4. Update parameters
        theta_new <- theta + delta

        # 5. Check convergence
        max_change <- max(abs(delta))
        if (verbose) {
            cat(sprintf(
                "Iter %d: Max Delta = %.6f, Objective ~ %.6f\n",
                iter, max_change, stats$gamma3^2 / stats$denom
            ))
        }

        if (max_change < tol) {
            theta <- theta_new
            if (verbose) cat("Converged.\n")
            break
        }

        theta <- theta_new
    }

    # Final recomputation
    final_res <- fn_residuals(theta)
    final_stats <- compute_pmm2_components(final_res)

    list(
        coefficients = theta,
        residuals = final_res,
        iterations = iter,
        converged = (iter < max_iter),
        pmm_stats = final_stats,
        method = "Iterative PMM2"
    )
}


#' Universal PMM2 estimator (One-step Global)
#'
#' Applies a single PMM2 correction to classical estimation results.
#'
#' @param theta_classical Parameter estimates from a classical method (e.g., MLE)
#' @param fn_residuals Function(theta) returning the residual vector
#' @param fn_jacobian Function(theta) returning the Jacobian matrix.
#'                    If NULL, numerical Jacobian via numDeriv is used
#' @param verbose Print progress
#'
#' @return List with results (theta, residuals, etc.)
#' @export
pmm2_nonlinear_onestep <- function(theta_classical, fn_residuals, fn_jacobian = NULL, verbose = FALSE) {
    if (verbose) cat("Starting One-step PMM2 correction...\n")

    # 1. Compute residuals and Jacobian at the classical estimate
    res <- fn_residuals(theta_classical)

    if (is.null(fn_jacobian)) {
        J <- compute_numerical_jacobian(fn_residuals, theta_classical)
    } else {
        J <- fn_jacobian(theta_classical)
    }

    # 2. PMM2 statistics
    stats <- compute_pmm2_components(res)
    if (is.null(stats)) {
        warning("PMM2 stats computation failed.")
        return(list(coefficients = theta_classical, method = "One-step PMM2 (Failed)"))
    }

    # 3. Compute single update step
    delta <- solve_pmm2_step(res, J, stats)

    # 4. Update
    theta_new <- theta_classical + delta

    if (verbose) {
        cat("Classical Theta:", paste(round(theta_classical, 4), collapse = ", "), "\n")
        cat("PMM2 Correction:", paste(round(delta, 4), collapse = ", "), "\n")
        cat("Final Theta:    ", paste(round(theta_new, 4), collapse = ", "), "\n")
    }

    # Final residuals
    final_res <- fn_residuals(theta_new)

    list(
        coefficients = theta_new,
        residuals = final_res,
        correction = delta,
        pmm_stats = stats,
        method = "One-step PMM2"
    )
}
