# pmm2_common.R - Common utilities for all PMM2 models


#' Universal PMM2 algorithm for all model types
#'
#' @param b_init Initial parameter estimates
#' @param X Design matrix
#' @param y Response vector
#' @param m2,m3,m4 Central moments
#' @param max_iter Maximum number of iterations
#' @param tol Tolerance for convergence
#' @param regularize Whether to add regularization
#' @param reg_lambda Regularization parameter
#' @param verbose Whether to print progress information
#' @param poly_terms Pre-computed polynomial coefficients (list with elements \code{A}, \code{B}, \code{C});
#'   allows passing custom values for special scenarios, otherwise they are computed from moments
#'
#' @return List with estimation results
#' @keywords internal
pmm2_algorithm <- function(b_init, X, y, m2, m3, m4,
                           max_iter = 50, tol = 1e-6,
                           regularize = TRUE, reg_lambda = 1e-8,
                           verbose = FALSE,
                           poly_terms = NULL) {
  # Current parameter estimates
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Calculate PMM2 polynomial coefficients
  # Derivation corresponds to equation (10) from the paper:
  # A = c3 * sigma^3, B = (c4 + 3) * sigma^4 - sigma^4 - 2 c3 * sigma^3 * y,
  # C = c3 * y^2 * sigma^3 - ((c4 + 3) * sigma^4 - sigma^4) * y - sigma^2 * c3 * sigma^3.
  # Substituting sigma^2 = m2, m3 = c3 sigma^3, m4 = (c4 + 3) sigma^4 gives the formulas below.
  if (is.null(poly_terms)) {
    A <- m3
    B <- m4 - m2^2 - 2 * m3 * y
    C <- m3 * y^2 - (m4 - m2^2) * y - m2 * m3
  } else {
    A <- poly_terms$A
    B <- poly_terms$B
    C <- poly_terms$C
    if (length(B) != nrow(X) || length(C) != nrow(X)) {
      stop("poly_terms$B and poly_terms$C must have length equal to the number of rows in X")
    }
  }

  # Track convergence history if verbose
  if (verbose) {
    conv_history <- numeric(max_iter)
  }

  # Main iteration loop
  for (iter in seq_len(max_iter)) {
    iterations <- iter

    # Calculate predicted values
    y_pred <- as.vector(X %*% b_cur)

    # Calculate Z1 = A*y_pred^2 + B*y_pred + C
    Z1 <- A*(y_pred^2) + B*y_pred + C

    # Form Z vector for each parameter
    p <- length(b_cur)
    Z <- numeric(p)
    for (r in 1:p) {
      Z[r] <- sum(Z1 * X[, r])
    }

    # Calculate derivative JZ11 = 2*A*y_pred + B
    JZ11 <- 2 * A * y_pred + B

    # Form Jacobian matrix
    JZs <- matrix(0, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        JZs[i, j] <- sum(JZ11 * X[, i] * X[, j])
      }
    }

    # Add regularization if needed
    if (regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Solve system JZs * delta = Z
    delta <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if (verbose) {
        cat("Error solving linear system:", conditionMessage(e), "\n")
        cat("Adding stronger regularization\n")
      }
      diag(JZs) <- diag(JZs) + 1e-4
      solve(JZs, Z)
    })

    # Update parameters
    b_new <- b_cur - delta
    diff_par <- sqrt(sum((b_new - b_cur)^2))

    # Store convergence history if verbose
    if (verbose) {
      conv_history[iter] <- diff_par
      if (iter %% 5 == 0 || iter == 1) {
        cat("Iteration", iter, ": Parameter change =",
            formatC(diff_par, digits = 8), "\n")
      }
    }

    b_cur <- b_new

    # Check convergence
    if (diff_par < tol) {
      converged <- TRUE
      if (verbose) cat("Converged after", iter, "iterations\n")
      break
    }
  }

  # Warning if maximum iterations reached without convergence
  if (!converged && verbose) {
    cat("Warning: Algorithm did not converge after", max_iter, "iterations\n")
  }

  # Calculate final residuals
  final_res <- as.numeric(y - X %*% b_cur)

  # Plot convergence history if verbose
  if (verbose && iterations > 1) {
    if (requireNamespace("graphics", quietly = TRUE)) {
      graphics::plot(1:iterations, conv_history[1:iterations], type = "b",
                     xlab = "Iteration", ylab = "Parameter change",
                     main = "Convergence history")
      graphics::abline(h = tol, col = "red", lty = 2)
    }
  }

  # Return results
  list(
    b = as.numeric(b_cur),
    convergence = converged,
    iterations = iterations,
    residuals = final_res
  )
}
