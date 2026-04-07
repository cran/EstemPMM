#' Compute MA innovations (forward recursion)
#'
#' @param x Time series (centered)
#' @param theta MA coefficients
#' @param q MA order
#' @return Vector of innovations
#' @keywords internal
ma_compute_innovations <- function(x, theta, q) {
    n <- length(x)
    innovations <- numeric(n)
    history <- rep(0, q)

    for (t in seq_len(n)) {
        ma_component <- if (q > 0) sum(theta * history) else 0
        innovations[t] <- x[t] - ma_component

        if (q > 0) {
            history <- c(innovations[t], history)[seq_len(q)]
        }
    }

    innovations
}


#' Compute SMA innovations (forward recursion)
#'
#' @param x Time series (centered)
#' @param Theta SMA coefficients
#' @param Q SMA order
#' @param s Seasonal period
#' @return Vector of innovations
#' @keywords internal
sma_compute_innovations <- function(x, Theta, Q, s) {
    n <- length(x)
    innovations <- numeric(n)

    for (t in seq_len(n)) {
        sma_component <- 0
        if (Q > 0) {
            for (J in 1:Q) {
                lag <- J * s
                if (t - lag >= 1) {
                    sma_component <- sma_component + Theta[J] * innovations[t - lag]
                }
            }
        }
        innovations[t] <- x[t] - sma_component
    }

    innovations
}


#' Build design matrix for MA
#'
#' @param intercept Intercept (mean)
#' @param residuals Initial residuals from CSS
#' @param x Time series
#' @param q MA order
#' @return List with X (design matrix) and y (response)
#' @keywords internal
ma_build_design <- function(intercept, residuals, x, q) {
    idx <- seq.int(q + 1L, length(x))
    X <- matrix(1, nrow = length(idx), ncol = q + 1L)

    for (j in seq_len(q)) {
        X[, j + 1L] <- residuals[idx - j]
    }

    y <- x[idx] - intercept

    list(X = X, y = y)
}


#' Build design matrix for SMA
#'
#' @param intercept Intercept (mean)
#' @param residuals Initial residuals from CSS
#' @param x Time series
#' @param Q SMA order
#' @param s Seasonal period
#' @return List with X (design matrix) and y (response)
#' @keywords internal
sma_build_design <- function(intercept, residuals, x, Q, s) {
    max_lag <- Q * s
    idx <- seq.int(max_lag + 1L, length(x))
    X <- matrix(1, nrow = length(idx), ncol = Q + 1L)

    for (J in seq_len(Q)) {
        lag_seasonal <- J * s
        X[, J + 1L] <- residuals[idx - lag_seasonal]
    }

    y <- x[idx] - intercept

    list(X = X, y = y)
}


#' Compute combined MA+SMA innovations (forward recursion)
#'
#' @param x Time series (centered)
#' @param theta_ma MA coefficients (length q)
#' @param theta_sma SMA coefficients (length Q)
#' @param q MA order
#' @param Q SMA order
#' @param s Seasonal period
#' @param multiplicative Logical, if TRUE use multiplicative SARIMA model
#' @return Vector of innovations
#' @keywords internal
ma_sma_compute_innovations <- function(x, theta_ma, theta_sma, q, Q, s, multiplicative = TRUE) {
    n <- length(x)
    innovations <- numeric(n)

    # History buffers
    ma_history <- if (q > 0) rep(0, q) else NULL

    for (t in seq_len(n)) {
        # Regular MA component
        ma_component <- if (q > 0) sum(theta_ma * ma_history) else 0

        # Seasonal MA component
        sma_component <- 0
        if (Q > 0) {
            for (J in 1:Q) {
                lag <- J * s
                if (t - lag >= 1) {
                    sma_component <- sma_component + theta_sma[J] * innovations[t - lag]
                }
            }
        }

        # Interaction component (Multiplicative only)
        interaction_component <- 0
        if (multiplicative && q > 0 && Q > 0) {
            for (j in 1:q) {
                for (J in 1:Q) {
                    lag_interaction <- j + J * s
                    if (t - lag_interaction >= 1) {
                        interaction_component <- interaction_component +
                            (theta_ma[j] * theta_sma[J]) * innovations[t - lag_interaction]
                    }
                }
            }
        }

        # Total innovation
        innovations[t] <- x[t] - ma_component - sma_component - interaction_component

        # Update MA history
        if (q > 0) {
            ma_history <- c(innovations[t], ma_history)[seq_len(q)]
        }
    }

    innovations
}


#' Build design matrix for combined MA+SMA model
#'
#' @param intercept Intercept (mean)
#' @param residuals Initial residuals from CSS
#' @param x Time series
#' @param q MA order
#' @param Q SMA order
#' @param s Seasonal period
#' @param multiplicative Logical, if TRUE use multiplicative SARIMA model
#' @return List with X (design matrix) and y (response)
#' @keywords internal
ma_sma_build_design <- function(intercept, residuals, x, q, Q, s, multiplicative = TRUE) {
    # Maximum lag determines where we start
    max_ma_lag <- q
    max_sma_lag <- Q * s
    max_lag <- max(max_ma_lag, max_sma_lag)

    if (multiplicative && q > 0 && Q > 0) {
        max_lag <- max(max_lag, q + Q * s)
    }

    idx <- seq.int(max_lag + 1L, length(x))

    # Design matrix: [intercept, ma_1, ..., ma_q, sma_1, ..., sma_Q, interaction_1_1, ...]
    num_interaction <- if (multiplicative && q > 0 && Q > 0) q * Q else 0
    ncol_total <- 1 + q + Q + num_interaction

    X <- matrix(1, nrow = length(idx), ncol = ncol_total)

    # Add regular MA lags
    if (q > 0) {
        for (j in seq_len(q)) {
            X[, j + 1L] <- residuals[idx - j]
        }
    }

    # Add seasonal MA lags
    if (Q > 0) {
        for (J in seq_len(Q)) {
            lag_seasonal <- J * s
            X[, q + J + 1L] <- residuals[idx - lag_seasonal]
        }
    }

    # Add interaction lags
    if (num_interaction > 0) {
        col_idx <- q + Q + 1L
        for (j in seq_len(q)) {
            for (J in seq_len(Q)) {
                lag_interaction <- j + J * s
                X[, col_idx + 1L] <- residuals[idx - lag_interaction]
                col_idx <- col_idx + 1L
            }
        }
    }

    y <- x[idx] - intercept

    list(X = X, y = y)
}


#' PMM2 solver (EstemPMM formula)
#'
#' @param b_init Initial coefficients vector (intercept and MA coefficients)
#' @param X Design matrix
#' @param Y Response vector
#' @param m2 Second central moment
#' @param m3 Third central moment
#' @param m4 Fourth central moment
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param verbose Print diagnostics
#'
#' @return List with coefficients, convergence, iterations
#' @keywords internal
ma_solve_pmm2 <- function(b_init, X, Y, m2, m3, m4,
                          max_iter = 50, tol = 1e-6,
                          verbose = FALSE) {
    b <- as.numeric(b_init)
    iterations <- 0L
    converged <- FALSE

    for (iter in seq_len(max_iter)) {
        iterations <- iter

        # S = X %*% b (predicted values)
        S <- as.vector(X %*% b)

        # Gradient Z1 (EstemPMM formula)
        Z1 <- m3 * S^2 +
            (m4 - m2^2 - 2 * m3 * Y) * S +
            (m3 * Y^2 - (m4 - m2^2) * Y - m2 * m3)

        # Z = t(X) %*% Z1
        Z <- as.numeric(t(X) %*% Z1)

        # Jacobian JZ11
        JZ11 <- 2 * m3 * S + (m4 - m2^2 - 2 * m3 * Y)

        # J = t(X) %*% (X * JZ11)
        J <- t(X) %*% (X * JZ11)

        # Newton step
        step <- tryCatch(solve(J, Z), error = function(e) NULL)

        if (is.null(step)) {
            if (verbose) message("  System is singular at iteration ", iter)
            break
        }

        # Update
        b_new <- b - step

        # Check convergence
        if (sqrt(sum((b_new - b)^2)) < tol) {
            b <- b_new
            converged <- TRUE
            if (verbose) message("  Converged at iteration ", iter)
            break
        }

        b <- b_new
    }

    list(
        coefficients = b[-1], # Remove intercept
        intercept = b[1],
        convergence = converged,
        iterations = iterations
    )
}


#' EstemPMM-style PMM2 estimator for MA models
#'
#' @param x Time series
#' @param q MA order
#' @param include.mean Include intercept
#' @param max_iter Maximum PMM2 iterations
#' @param verbose Print diagnostics
#'
#' @return List with ma_coef, mean, innovations, convergence, method
#' @keywords internal
estpmm_style_ma <- function(x,
                            q = 1,
                            include.mean = TRUE,
                            max_iter = 50,
                            verbose = FALSE) {
    n <- length(x)

    # STEP 1: CSS fit for initial estimates
    css_fit <- tryCatch(
        {
            stats::arima(x,
                order = c(0, 0, q),
                method = "CSS",
                include.mean = include.mean
            )
        },
        error = function(e) {
            stats::arima(x,
                order = c(0, 0, q),
                method = "CSS-ML",
                include.mean = include.mean
            )
        }
    )

    # Extract initial parameters
    coefs_css <- coef(css_fit)

    theta_init <- numeric(q)
    for (j in 1:q) {
        coef_name <- paste0("ma", j)
        if (coef_name %in% names(coefs_css)) {
            theta_init[j] <- coefs_css[coef_name]
        }
    }

    intercept_init <- if (include.mean && "intercept" %in% names(coefs_css)) {
        coefs_css["intercept"]
    } else {
        0
    }

    # Compute residuals from CSS
    residuals_css <- as.numeric(residuals(css_fit))

    # STEP 2: Build design matrix from FIXED residuals
    design <- ma_build_design(intercept_init, residuals_css, x, q)
    X <- design$X
    y <- design$y

    # STEP 3: Compute moments from CSS residuals
    # Use effective residuals (after removing initial q observations)
    eff_residuals <- residuals_css[(q + 1):n]

    m2 <- mean(eff_residuals^2)
    m3 <- mean(eff_residuals^3)
    m4 <- mean(eff_residuals^4)

    # STEP 4: PMM2 optimization
    b_init <- c(intercept_init, theta_init)

    pmm2_result <- ma_solve_pmm2(b_init, X, y, m2, m3, m4,
        max_iter = max_iter,
        tol = 1e-6,
        verbose = verbose
    )

    # Final parameters
    theta_final <- pmm2_result$coefficients
    intercept_final <- pmm2_result$intercept

    # STEP 5: Compute final innovations
    x_centered <- as.numeric(x) - intercept_final
    innovations <- ma_compute_innovations(x_centered, theta_final, q)

    list(
        ma_coef = theta_final,
        mean = intercept_final,
        innovations = innovations,
        convergence = pmm2_result$convergence,
        iterations = pmm2_result$iterations,
        css_estimates = theta_init,
        method = "EstemPMM-style PMM2"
    )
}


#' EstemPMM-style PMM2 estimator for SMA models
#'
#' @param x Time series
#' @param Q SMA order
#' @param s Seasonal period
#' @param include.mean Include intercept
#' @param max_iter Maximum PMM2 iterations
#' @param verbose Print diagnostics
#'
#' @return List with sma_coef, mean, innovations, convergence, method
#' @keywords internal
estpmm_style_sma <- function(x,
                             Q = 1,
                             s = 4,
                             include.mean = TRUE,
                             max_iter = 50,
                             verbose = FALSE) {
    n <- length(x)

    # STEP 1: CSS fit
    css_fit <- tryCatch(
        {
            stats::arima(x,
                order = c(0, 0, 0),
                seasonal = list(order = c(0, 0, Q), period = s),
                method = "CSS",
                include.mean = include.mean
            )
        },
        error = function(e) {
            stats::arima(x,
                order = c(0, 0, 0),
                seasonal = list(order = c(0, 0, Q), period = s),
                method = "CSS-ML",
                include.mean = include.mean
            )
        }
    )

    coefs_css <- coef(css_fit)

    Theta_init <- numeric(Q)
    for (J in 1:Q) {
        coef_name <- paste0("sma", J)
        if (coef_name %in% names(coefs_css)) {
            Theta_init[J] <- coefs_css[coef_name]
        }
    }

    intercept_init <- if (include.mean && "intercept" %in% names(coefs_css)) {
        coefs_css["intercept"]
    } else {
        0
    }

    residuals_css <- as.numeric(residuals(css_fit))

    # STEP 2: Build design matrix
    design <- sma_build_design(intercept_init, residuals_css, x, Q, s)
    X <- design$X
    y <- design$y

    # STEP 3: Compute moments
    max_lag <- Q * s
    eff_residuals <- residuals_css[(max_lag + 1):n]

    m2 <- mean(eff_residuals^2)
    m3 <- mean(eff_residuals^3)
    m4 <- mean(eff_residuals^4)

    # STEP 4: PMM2 optimization
    b_init <- c(intercept_init, Theta_init)

    pmm2_result <- ma_solve_pmm2(b_init, X, y, m2, m3, m4,
        max_iter = max_iter,
        tol = 1e-6,
        verbose = verbose
    )

    Theta_final <- pmm2_result$coefficients
    intercept_final <- pmm2_result$intercept

    # STEP 5: Compute final innovations
    x_centered <- as.numeric(x) - intercept_final
    innovations <- sma_compute_innovations(x_centered, Theta_final, Q, s)

    list(
        sma_coef = Theta_final,
        mean = intercept_final,
        innovations = innovations,
        convergence = pmm2_result$convergence,
        iterations = pmm2_result$iterations,
        css_estimates = Theta_init,
        method = "EstemPMM-style PMM2"
    )
}


#' EstemPMM-style PMM2 estimator for combined MA+SMA models
#'
#' @param x Time series
#' @param q MA order
#' @param Q SMA order
#' @param s Seasonal period
#' @param include.mean Include intercept
#' @param max_iter Maximum PMM2 iterations
#' @param verbose Print diagnostics
#' @param multiplicative Logical, if TRUE use multiplicative SARIMA model
#'
#' @return List with ma_coef, sma_coef, mean, innovations, convergence, method
#' @keywords internal
estpmm_style_ma_sma <- function(x,
                                q = 1,
                                Q = 1,
                                s = 12,
                                include.mean = TRUE,
                                max_iter = 50,
                                verbose = FALSE,
                                multiplicative = TRUE) {
    n <- length(x)

    # Validate inputs
    if (q < 1) stop("q must be at least 1")
    if (Q < 1) stop("Q must be at least 1")
    if (s < 2) stop("Seasonal period must be at least 2")

    # STEP 1: CSS fit for initial estimates
    css_fit <- tryCatch(
        {
            stats::arima(x,
                order = c(0, 0, q),
                seasonal = list(order = c(0, 0, Q), period = s),
                method = "CSS",
                include.mean = include.mean
            )
        },
        error = function(e) {
            stats::arima(x,
                order = c(0, 0, q),
                seasonal = list(order = c(0, 0, Q), period = s),
                method = "CSS-ML",
                include.mean = include.mean
            )
        }
    )

    # Extract initial parameters
    coefs_css <- coef(css_fit)

    # Extract MA coefficients
    theta_ma_init <- numeric(q)
    for (j in 1:q) {
        coef_name <- paste0("ma", j)
        if (coef_name %in% names(coefs_css)) {
            theta_ma_init[j] <- coefs_css[coef_name]
        }
    }

    # Extract SMA coefficients
    theta_sma_init <- numeric(Q)
    for (J in 1:Q) {
        coef_name <- paste0("sma", J)
        if (coef_name %in% names(coefs_css)) {
            theta_sma_init[J] <- coefs_css[coef_name]
        }
    }

    intercept_init <- if (include.mean && "intercept" %in% names(coefs_css)) {
        coefs_css["intercept"]
    } else {
        0
    }

    # Compute residuals from CSS
    residuals_css <- as.numeric(residuals(css_fit))

    # STEP 2: Build design matrix from FIXED residuals
    design <- ma_sma_build_design(intercept_init, residuals_css, x, q, Q, s, multiplicative)
    X <- design$X
    y <- design$y

    # STEP 3: Compute moments from CSS residuals
    max_lag <- max(q, Q * s)
    if (multiplicative) {
        max_lag <- max(max_lag, q + Q * s)
    }
    eff_residuals <- residuals_css[(max_lag + 1):n]

    m2 <- mean(eff_residuals^2)
    m3 <- mean(eff_residuals^3)
    m4 <- mean(eff_residuals^4)

    # STEP 4: PMM2 optimization
    # Combined parameter vector: [intercept, ma_1, ..., ma_q, sma_1, ..., sma_Q, interaction_1_1, ...]
    b_init <- c(intercept_init, theta_ma_init, theta_sma_init)

    # Add initial interaction terms if multiplicative
    if (multiplicative) {
        interaction_init <- numeric(q * Q)
        idx <- 1
        for (j in 1:q) {
            for (J in 1:Q) {
                interaction_init[idx] <- theta_ma_init[j] * theta_sma_init[J]
                idx <- idx + 1
            }
        }
        b_init <- c(b_init, interaction_init)
    }

    pmm2_result <- ma_solve_pmm2(b_init, X, y, m2, m3, m4,
        max_iter = max_iter,
        tol = 1e-6,
        verbose = verbose
    )

    # Extract final parameters
    all_coefs <- pmm2_result$coefficients
    theta_ma_final <- if (q > 0) all_coefs[1:q] else numeric(0)
    theta_sma_final <- if (Q > 0) all_coefs[(q + 1):(q + Q)] else numeric(0)

    # We ignore the estimated interaction coefficients for the final result,
    # relying on the multiplicative structure implied by theta_ma and theta_sma

    intercept_final <- pmm2_result$intercept

    # STEP 5: Compute final innovations
    x_centered <- as.numeric(x) - intercept_final
    innovations <- ma_sma_compute_innovations(x_centered, theta_ma_final, theta_sma_final, q, Q, s, multiplicative)

    list(
        ma_coef = theta_ma_final,
        sma_coef = theta_sma_final,
        mean = intercept_final,
        innovations = innovations,
        convergence = pmm2_result$convergence,
        iterations = pmm2_result$iterations,
        css_estimates = list(ma = theta_ma_init, sma = theta_sma_init),
        method = "EstemPMM-style PMM2"
    )
}
