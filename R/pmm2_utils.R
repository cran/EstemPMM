# pmm2_utils.R

#' Universal solver for PMM2 system of equations
#'
#' @param b_init Initial parameter estimates (usually from OLS or MLE)
#' @param X Design matrix (including intercept and all predictors)
#' @param y Response vector
#' @param m2 Second central moment of residuals
#' @param m3 Third central moment of residuals
#' @param m4 Fourth central moment of residuals
#' @param max_iter Maximum number of iterations
#' @param tol Convergence tolerance
#' @param regularize Whether to add regularization to Jacobian matrix
#' @param reg_lambda Regularization parameter
#' @param verbose Print progress information
#'
#' @return Vector of PMM2 parameter estimates
#' @keywords internal
solve_pmm2 <- function(b_init, X, y, m2, m3, m4,
                       max_iter = 1000, tol = 1e-5,
                       regularize = TRUE, reg_lambda = 1e-8,
                       verbose = FALSE) {

  # Get dimensions
  P <- length(b_init)
  n <- nrow(X)

  # Define coefficients for PMM2 polynomial
  A <- m3
  B <- m4 - m2^2 - 2*m3*y
  C <- m3*y^2 - (m4 - m2^2)*y - m2*m3

  # Current parameter estimates
  b_pmm2 <- b_init

  # Initialize Jacobian matrix
  JZs <- matrix(0, nrow = P, ncol = P)

  # Iteration loop
  for (iter in 1:max_iter) {
    # Calculate predicted values
    y_pred <- as.numeric(X %*% b_pmm2)

    # Calculate Z1 = A*y_pred^2 + B*y_pred + C
    Z1 <- A*y_pred^2 + B*y_pred + C

    # Form Z vector for each parameter
    Z <- numeric(P)
    for (r in 1:P) {
      Z[r] <- sum(Z1 * X[, r])
    }

    # Calculate derivative JZ11 = 2*A*y_pred + B
    JZ11 <- 2*A*y_pred + B

    # Form Jacobian matrix
    for (i in 1:P) {
      for (j in 1:P) {
        JZs[i, j] <- sum(JZ11 * X[, i] * X[, j])
      }
    }

    # Add regularization if requested
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
    b_pmm2 <- b_pmm2 - delta

    # Check convergence
    if (sqrt(sum(delta^2)) < tol) {
      if (verbose) {
        cat("Converged after", iter, "iterations\n")
      }
      break
    }

    # Print progress if requested
    if (verbose && iter %% 10 == 0) {
      cat("Iteration", iter, ", change = ", sqrt(sum(delta^2)), "\n")
    }
  }

  # Warning if not converged
  if (iter == max_iter && verbose) {
    cat("Warning: Maximum iterations reached without convergence\n")
  }

  return(b_pmm2)
}


#' PMM2 fitting algorithm - unified implementation
#'
#' Core iterative algorithm for PMM2 parameter estimation
#'
#' @param b_init initial parameter estimates (typically from OLS)
#' @param X design matrix
#' @param y response vector
#' @param m2,m3,m4 central moments
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param regularize logical, add small value to diagonal for numerical stability
#' @param reg_lambda regularization parameter (if regularize=TRUE)
#' @param verbose logical, whether to print progress information
#'
#' @return A list with components:
#'   \item{b}{estimated parameters}
#'   \item{convergence}{logical convergence status}
#'   \item{iterations}{number of iterations performed}
#'
#' @keywords internal
.pmm2_fit <- function(b_init, X, y,
                      m2, m3, m4,
                      max_iter=50, tol=1e-6,
                      regularize=TRUE, reg_lambda=1e-8,
                      verbose=FALSE) {
  n <- nrow(X)
  p <- ncol(X)

  # Compute A, B, C for the entire dataset at once
  A <- m3
  B <- (m4 - m2^2) - 2*m3*y
  C <- m3*(y^2) - ((m4 - m2^2)*y) - m2*m3

  # Current parameter estimates
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Track convergence history if verbose
  if(verbose) {
    conv_history <- numeric(max_iter)
  }

  for (iter in seq_len(max_iter)) {
    iterations <- iter

    # Compute Yx = X %*% b_cur (predictor)
    # More efficient matrix multiplication
    Yx <- as.vector(X %*% b_cur)  # Ensure Yx is a vector

    # Compute Z1 = A*Yx^2 + B*Yx + C
    Z1 <- A*(Yx^2) + B*Yx + C

    # Form Z - vector of equations
    # Improved handling for multi-variable cases
    if (p > 1) {
      # Perform crossprod for all predictors except intercept
      Z_rest <- crossprod(X[, -1, drop=FALSE], Z1)

      # Check if Z_rest is a matrix (more than 1 predictor)
      if (is.matrix(Z_rest)) {
        Z_rest <- as.vector(Z_rest)
      }

      # Combine results
      Z <- c(sum(Z1), Z_rest)
    } else {
      # If only intercept, Z is just the sum
      Z <- sum(Z1)
    }

    # Form JZs - Jacobian matrix
    JZ11 <- 2*A*Yx + B

    # Vectorized computation of Jacobian matrix
    JZs <- matrix(NA, p, p)
    JZs[1,1] <- sum(JZ11)

    # First row and first column of Jacobian
    for (ii in 2:p) {
      tmp <- JZ11 * X[,ii]
      JZs[1,ii] <- sum(tmp)
      JZs[ii,1] <- JZs[1,ii]  # Jacobian is symmetric
    }

    # Remaining elements of Jacobian
    for (ii in 2:p) {
      for (jj in 2:p) {
        tmp <- JZ11 * X[,ii] * X[,jj]
        JZs[ii,jj] <- sum(tmp)
      }
    }

    # Apply regularization if needed to avoid singularity
    if (regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Solve system JZs * delta = Z
    step <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if(verbose) {
        cat("Error solving linear system in iteration", iter, ":", conditionMessage(e), "\n")
      }
      # Fallback to pseudoinverse or more stable solver
      if(requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(JZs) %*% Z
      } else {
        warning("Failed to solve linear system. Consider installing 'MASS' package.")
        rep(NA, length(Z))
      }
    })

    # Check for numerical problems
    if (any(is.na(step)) || any(is.infinite(step))) {
      warning("Numerical problems encountered in iteration ", iter,
              ". Algorithm stopped.")
      break
    }

    # Update parameters
    b_new <- b_cur - step
    diff_par <- sqrt(sum((b_new - b_cur)^2))

    # Store convergence history if verbose
    if(verbose) {
      conv_history[iter] <- diff_par
      if(iter %% 5 == 0 || iter == 1) {
        cat("Iteration", iter, ": Parameter change =",
            formatC(diff_par, digits=8), "\n")
      }
    }

    b_cur <- b_new

    # Check convergence
    if (diff_par < tol) {
      converged <- TRUE
      if(verbose) cat("Converged after", iter, "iterations\n")
      break
    }
  }

  # Warning if max iterations reached without convergence
  if(!converged && verbose) {
    cat("Warning: Algorithm did not converge after", max_iter, "iterations\n")
  }

  # Plot convergence history if verbose
  if(verbose && iterations > 1) {
    if(requireNamespace("graphics", quietly = TRUE)) {
      graphics::plot(1:iterations, conv_history[1:iterations], type="b",
                     xlab="Iteration", ylab="Parameter change",
                     main="Convergence history")
      graphics::abline(h=tol, col="red", lty=2)
    }
  }

  # Return results - ensure b is a numeric vector, not a matrix
  list(
    b = as.numeric(b_cur),
    convergence = converged,
    iterations = iterations
  )
}

#' PMM2 fitting algorithm for time series models
#'
#' Implementation of PMM2 algorithm for time series models with various approaches
#' depending on the model structure (AR, MA, ARMA, ARIMA)
#'
#' @param b_init initial parameter estimates (AR followed by MA)
#' @param innovations_init initial innovations (errors)
#' @param model_info list with model structure information
#' @param m2,m3,m4 central moments
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param regularize logical, add small value to diagonal for numerical stability
#' @param reg_lambda regularization parameter (if regularize=TRUE)
#' @param verbose logical, whether to print progress information
#'
#' @return A list with components:
#'   \item{b}{estimated parameters}
#'   \item{convergence}{logical convergence status}
#'   \item{iterations}{number of iterations performed}
#'   \item{innovations}{final innovations (errors)}
#'
#' @keywords internal
.ts_pmm2_fit <- function(b_init, innovations_init, model_info,
                         m2, m3, m4,
                         max_iter = 50, tol = 1e-6,
                         regularize = TRUE, reg_lambda = 1e-8,
                         verbose = FALSE) {

  # Extract model information
  x <- model_info$x
  model_type <- model_info$model_type
  ar_order <- model_info$ar_order
  ma_order <- model_info$ma_order

  n <- length(x)
  total_params <- ar_order + ma_order

  # Current parameter estimates
  b_cur <- b_init
  converged <- FALSE
  iterations <- 0

  # Current innovations
  innovations <- innovations_init

  # Check for non-finite values in innovations
  if(any(is.infinite(innovations)) || any(is.na(innovations))) {
    warning("Non-finite values in initial innovations. Replacing with estimated values.")
    bad_idx <- is.infinite(innovations) | is.na(innovations)
    if(sum(!bad_idx) > 0) {
      innovations[bad_idx] <- mean(innovations[!bad_idx])
    } else {
      innovations[bad_idx] <- 0
    }
  }

  # Track convergence history if verbose
  if(verbose) {
    conv_history <- numeric(max_iter)
  }

  for(iter in seq_len(max_iter)) {
    iterations <- iter

    # Update innovations based on model type and current parameters
    if(model_type == "ar") {
      # For AR models, calculate predicted values
      X <- create_ar_matrix(x, ar_order)
      y <- x[(ar_order + 1):n]

      # Verify dimensions to prevent errors
      if(nrow(X) != length(y)) {
        warning("Dimension mismatch in AR model. Adjusting matrices.")
        min_len <- min(nrow(X), length(y))
        X <- X[1:min_len, , drop = FALSE]
        y <- y[1:min_len]
      }

      y_pred <- as.vector(X %*% b_cur)

      # Form matrices for the PMM algorithm
      A <- m3
      B <- (m4 - m2^2) - 2*m3*y
      C <- m3*(y^2) - ((m4 - m2^2)*y) - m2*m3

      # Calculate Z1 = A*yhat^2 + B*yhat + C
      Z1 <- A*(y_pred^2) + B*y_pred + C

      # Form Z vector for each parameter
      Z <- numeric(ar_order)
      for(i in 1:ar_order) {
        Z[i] <- sum(Z1 * X[, i])
      }

      # Form Jacobian matrix
      JZ11 <- 2*A*y_pred + B
      JZs <- matrix(0, ar_order, ar_order)

      for(i in 1:ar_order) {
        for(j in 1:ar_order) {
          JZs[i, j] <- sum(JZ11 * X[, i] * X[, j])
        }
      }

      # Update innovations
      innovations <- y - y_pred

    } else if(model_type %in% c("ma", "arma", "arima")) {
      # Extract AR and MA parts from current parameters
      ar_part <- if(ar_order > 0) b_cur[1:ar_order] else numeric(0)
      ma_part <- if(ma_order > 0) b_cur[(ar_order+1):total_params] else numeric(0)

      # For MA-based models, we need to update innovations using ARMA framework
      if(model_type == "ma") {
        # For pure MA, update innovations directly
        innovations <- tryCatch({
          update_ma_innovations(x, ma_part)
        }, error = function(e) {
          if(verbose) {
            cat("Error updating MA innovations:", conditionMessage(e), "\n")
          }
          # Return previous innovations if update fails
          innovations_init
        })
      } else {
        # For ARMA and ARIMA, use arima with fixed parameters
        arima_order <- c(ar_order, 0, ma_order) # d=0 as we're working with (differenced) x

        # Handle potential errors in model fitting
        curr_model <- tryCatch({
          stats::arima(x, order = arima_order,
                       include.mean = FALSE,
                       fixed = c(ar_part, ma_part))
        }, error = function(e) {
          if(verbose) {
            cat("Error in model fitting:", conditionMessage(e), "\n")
          }
          return(NULL)
        })

        if(is.null(curr_model)) {
          warning("Failed to update model in iteration ", iter, ". Using alternative approach.")
          # Use previous innovations with small adjustment
          innovations <- innovations * (1 + rnorm(length(innovations), 0, 0.01))
        } else {
          innovations <- residuals(curr_model)
        }
      }

      # Check for non-finite values in innovations
      if(any(is.infinite(innovations)) || any(is.na(innovations))) {
        warning("Non-finite innovations detected in iteration ", iter, ". Regularizing values.")
        bad_idx <- is.infinite(innovations) | is.na(innovations)
        if(sum(!bad_idx) > 0) {
          innovations[bad_idx] <- mean(innovations[!bad_idx])
        } else {
          innovations[bad_idx] <- 0
        }
      }

      # Prepare design matrices for both AR and MA parts
      J <- matrix(0, nrow = n, ncol = total_params)

      # For AR part, use lagged values of the series
      if(ar_order > 0) {
        for(i in 1:ar_order) {
          # Ensure indices are valid
          if(i+1 <= n && 1 <= n-i) {
            J[(i+1):n, i] <- x[1:(n-i)]
          }
        }
      }

      # For MA part, use lagged innovations
      if(ma_order > 0) {
        for(i in 1:ma_order) {
          # Ensure indices are valid
          if(i+1 <= n && 1 <= n-i && i <= length(innovations)) {
            J[(i+1):n, ar_order+i] <- innovations[1:(n-i)]
          }
        }
      }

      # Remove initial rows that can't be used for all lags
      start_idx <- max(ar_order, ma_order) + 1
      if(start_idx > 1 && start_idx <= n) {
        J_valid <- J[start_idx:n, , drop = FALSE]
        x_valid <- x[start_idx:n]
        innovations_valid <- innovations[start_idx:n]
      } else {
        # Handle the case when start_idx is too large
        warning("Not enough data points for the specified model orders. Using all available data.")
        J_valid <- J
        x_valid <- x
        innovations_valid <- innovations
      }

      # Ensure we have enough data
      if(nrow(J_valid) < 1) {
        warning("Insufficient data for PMM2 estimation. Stopping at iteration ", iter)
        break
      }

      # Compute PMM2 coefficients
      A <- m3
      B <- (m4 - m2^2) - 2*m3*x_valid
      C <- m3*(x_valid^2) - ((m4 - m2^2)*x_valid) - m2*m3

      # Form Z vector for each parameter
      Z <- numeric(total_params)
      for(j in 1:total_params) {
        Z1 <- A*(innovations_valid^2) + B*innovations_valid + C
        Z[j] <- sum(Z1 * J_valid[, j])
      }

      # Form Jacobian matrix
      JZ11 <- 2*A*innovations_valid + B
      JZs <- matrix(0, total_params, total_params)

      for(i in 1:total_params) {
        for(j in 1:total_params) {
          JZs[i, j] <- sum(JZ11 * J_valid[, i] * J_valid[, j])
        }
      }
    }

    # Apply regularization if needed
    if(regularize) {
      diag(JZs) <- diag(JZs) + reg_lambda
    }

    # Check if JZs is invertible
    if(any(is.na(JZs)) || any(is.infinite(JZs)) ||
       abs(det(JZs)) < 1e-10 || any(diag(JZs) < 1e-10)) {
      warning("Jacobian matrix is near-singular in iteration ", iter, ". Adding stronger regularization.")
      diag(JZs) <- diag(JZs) + 1e-4
    }

    # Solve system JZs * delta = Z
    step <- tryCatch({
      solve(JZs, Z)
    }, error = function(e) {
      if(verbose) {
        cat("Error solving linear system in iteration", iter, ":", conditionMessage(e), "\n")
      }
      # Fallback to pseudoinverse
      if(requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(JZs) %*% Z
      } else {
        warning("Failed to solve linear system. Consider installing 'MASS' package.")
        rep(NA, length(Z))
      }
    })

    # Check for numerical problems
    if(any(is.na(step)) || any(is.infinite(step))) {
      warning("Numerical problems encountered in iteration ", iter, ". Using scaled step.")
      # Replace NA/infinite values with small values
      step[is.na(step) | is.infinite(step)] <- 0.01 * sign(Z[is.na(step) | is.infinite(step)])
    }

    # Limit step size to prevent divergence
    max_step <- 0.5
    if(max(abs(step)) > max_step) {
      step <- step * (max_step / max(abs(step)))
      if(verbose) cat("Step size limited in iteration ", iter, "\n")
    }

    # Update parameters
    b_new <- b_cur - step
    diff_par <- sqrt(sum((b_new - b_cur)^2))

    # Store convergence history if verbose
    if(verbose) {
      conv_history[iter] <- diff_par
      if(iter %% 5 == 0 || iter == 1) {
        cat("Iteration", iter, ": Parameter change =",
            formatC(diff_par, digits = 8), "\n")
      }
    }

    b_cur <- b_new

    # Check convergence
    if(diff_par < tol) {
      converged <- TRUE
      if(verbose) cat("Converged after", iter, "iterations\n")
      break
    }
  }

  # Warning if max iterations reached without convergence
  if(!converged && verbose) {
    cat("Warning: Algorithm did not converge after", max_iter, "iterations\n")
  }

  # Final innovation update based on model type
  if(model_type == "ar") {
    X <- create_ar_matrix(x, ar_order)
    y <- x[(ar_order + 1):n]

    # Verify dimensions to prevent errors
    if(nrow(X) != length(y)) {
      min_len <- min(nrow(X), length(y))
      X <- X[1:min_len, , drop = FALSE]
      y <- y[1:min_len]
    }

    final_innovations <- y - X %*% b_cur
  } else {
    # Extract AR and MA parts
    ar_part <- if(ar_order > 0) b_cur[1:ar_order] else numeric(0)
    ma_part <- if(ma_order > 0) b_cur[(ar_order+1):total_params] else numeric(0)

    # Arima order for the final model
    arima_order <- c(ar_order, 0, ma_order)

    # Fit final model to get innovations
    final_model <- tryCatch({
      stats::arima(x, order = arima_order,
                   include.mean = FALSE,
                   fixed = c(ar_part, ma_part))
    }, error = function(e) {
      if(verbose) {
        cat("Error in final model fitting:", conditionMessage(e), "\n")
      }
      NULL
    })

    if(!is.null(final_model)) {
      final_innovations <- residuals(final_model)
    } else {
      final_innovations <- innovations
      warning("Failed to compute final innovations. Using last iteration values.")
    }
  }

  # Return results
  list(
    b = as.numeric(b_cur),
    convergence = converged,
    iterations = iterations,
    innovations = final_innovations
  )
}


#' Calculate kurtosis from data
#'
#' @param x numeric vector
#' @param excess logical, whether to return excess kurtosis (kurtosis - 3)
#'
#' @return Kurtosis value
#' @export
pmm_kurtosis <- function(x, excess = TRUE) {
  # Remove NAs
  x <- x[!is.na(x)]
  n <- length(x)

  if(n < 4) {
    warning("At least 4 non-missing values are needed to compute kurtosis")
    return(NA)
  }

  # Center the data
  x_centered <- x - mean(x)

  # Calculate moments
  m2 <- mean(x_centered^2)
  m4 <- mean(x_centered^4)

  # Calculate kurtosis
  kurt <- m4 / (m2^2)

  # Return excess kurtosis if requested
  if(excess) {
    kurt <- kurt - 3
  }

  return(kurt)
}

#' Calculate skewness from data
#'
#' @param x numeric vector
#'
#' @return Skewness value
#' @export
pmm_skewness <- function(x) {
  # Remove NAs
  x <- x[!is.na(x)]
  n <- length(x)

  if(n < 3) {
    warning("At least 3 non-missing values are needed to compute skewness")
    return(NA)
  }

  # Center the data
  x_centered <- x - mean(x)

  # Calculate moments
  m2 <- mean(x_centered^2)
  m3 <- mean(x_centered^3)

  # Calculate skewness
  skew <- m3 / (m2^(3/2))

  return(skew)
}


#' Calculate moments and cumulants of error distribution
#'
#' @param errors numeric vector of errors
#' @return list with moments, cumulants and theoretical variance reduction coefficient
#' @export
compute_moments <- function(errors) {
  m2 <- mean(errors^2)
  m3 <- mean(errors^3)
  m4 <- mean(errors^4)

  c3 <- m3 / m2^(3/2)  # Skewness coefficient
  c4 <- m4 / m2^2 - 3  # Excess kurtosis coefficient

  # Theoretical variance reduction coefficient
  g <- 1 - c3^2 / (2 + c4)

  return(list(m2 = m2, m3 = m3, m4 = m4,
              c3 = c3, c4 = c4,
              g = g))
}

#' Calculate theoretical skewness, kurtosis coefficients and variance reduction factor
#'
#' @param m2,m3,m4 central moments of second, third and fourth orders
#'
#' @return List with fields `c3`, `c4` and `g`
#' @export
pmm2_variance_factor <- function(m2, m3, m4) {
  if (is.na(m2) || m2 <= 0) {
    return(list(c3 = NA_real_, c4 = NA_real_, g = NA_real_))
  }
  c3 <- m3 / m2^(3/2)
  c4 <- m4 / m2^2 - 3
  denom <- 2 + c4
  g <- if (is.na(denom) || denom == 0) NA_real_ else 1 - (c3^2) / denom
  list(c3 = c3, c4 = c4, g = g)
}

#' Calculate theoretical variance matrices for OLS and PMM2
#'
#' @param X Design matrix with column of ones
#' @param m2,m3,m4 central moments of OLS residuals
#'
#' @return List with fields `ols`, `pmm2`, `c3`, `c4`, `g`
#' @export
pmm2_variance_matrices <- function(X, m2, m3, m4) {
  X <- as.matrix(X)
  XtX <- crossprod(X)
  V1 <- tryCatch({
    m2 * solve(XtX)
  }, error = function(e) {
    stop("Failed to invert matrix X'X: ", conditionMessage(e), call. = FALSE)
  })

  vf <- pmm2_variance_factor(m2, m3, m4)
  if (is.na(vf$g)) {
    V2 <- matrix(NA_real_, nrow = nrow(V1), ncol = ncol(V1))
  } else {
    V2 <- vf$g * V1
  }

  list(ols = V1, pmm2 = V2, c3 = vf$c3, c4 = vf$c4, g = vf$g)
}
