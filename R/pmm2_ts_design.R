# pmm2_ts_design.R - Functions for working with time series design matrices

#' Validate and prepare time series parameters
#'
#' @param x Time series data
#' @param order Model order specification
#' @param model_type Model type (ar, ma, arma, or arima)
#' @param include.mean Whether to include mean/intercept
#'
#' @return List of validated parameters and model information
#' @keywords internal
validate_ts_parameters <- function(x, order, model_type, include.mean) {
  # Check input data
  if (missing(x)) {
    stop("Missing argument 'x'")
  }

  # Convert input data to numeric vector
  x <- as.numeric(x)

  if (!is.numeric(x)) {
    stop("'x' must be a numeric vector")
  }

  # Check for NA and infinite values
  if (any(is.na(x)) || any(is.infinite(x))) {
    warning("NA or infinite values detected in input series. They will be removed.")
    x <- x[!is.na(x) & !is.infinite(x)]
    if (length(x) < 10) {
      stop("Too few valid observations after removing NA/infinite values")
    }
  }

  if (missing(order)) {
    stop("Missing argument 'order'")
  }

  # Parse order parameter depending on model_type
  if (model_type == "ar") {
    if (!is.numeric(order) || length(order) != 1)
      stop("For AR models 'order' must be a single integer")
    ar_order <- as.integer(order)
    ma_order <- 0
    d <- 0
    if (ar_order <= 0)
      stop("AR order must be positive")

    # Check if there is enough data
    if (length(x) <= ar_order + 1) {
      stop("Too few observations for AR model of order ", ar_order)
    }
  } else if (model_type == "ma") {
    if (!is.numeric(order) || length(order) != 1)
      stop("For MA models 'order' must be a single integer")
    ar_order <- 0
    ma_order <- as.integer(order)
    d <- 0
    if (ma_order <= 0)
      stop("MA order must be positive")

    # Check if there is enough data
    if (length(x) <= ma_order + 1) {
      stop("Too few observations for MA model of order ", ma_order)
    }
  } else if (model_type == "arma") {
    if (!is.numeric(order) || length(order) != 2)
      stop("For ARMA models 'order' must be a vector of length 2 (AR order, MA order)")
    ar_order <- as.integer(order[1])
    ma_order <- as.integer(order[2])
    d <- 0
    if (ar_order < 0 || ma_order < 0)
      stop("AR and MA orders must be non-negative")
    if (ar_order == 0 && ma_order == 0)
      stop("At least one of AR or MA orders must be positive")

    # Check if there is enough data
    if (length(x) <= max(ar_order, ma_order) + 1) {
      stop("Too few observations for ARMA model of orders (", ar_order, ",", ma_order, ")")
    }
  } else if (model_type == "arima") {
    if (!is.numeric(order) || length(order) != 3)
      stop("For ARIMA models 'order' must be a vector of length 3 (AR order, differencing, MA order)")
    ar_order <- as.integer(order[1])
    d <- as.integer(order[2])
    ma_order <- as.integer(order[3])
    if (ar_order < 0 || ma_order < 0 || d < 0)
      stop("AR, differencing, and MA orders must be non-negative")
    if (ar_order == 0 && ma_order == 0 && d == 0)
      stop("At least one of AR, differencing, or MA orders must be positive")

    # Check if there is enough data after differencing
    if (length(x) <= d + max(ar_order, ma_order) + 1) {
      stop("Too few observations for ARIMA model after differencing")
    }
  } else {
    stop("Unknown model type: ", model_type)
  }

  # Save original series
  orig_x <- as.numeric(x)

  list(
    original_x = orig_x,
    ar_order   = ar_order,
    ma_order   = ma_order,
    d          = d,
    model_type = model_type,
    include.mean = include.mean
  )
}

#' Create design matrix for AR model
#'
#' @param x Centered time series
#' @param p AR order
#' @return Design matrix with lagged values
#' @keywords internal
create_ar_matrix <- function(x, p) {
  n <- length(x)
  if (n <= p) {
    stop("Insufficient data points for AR order p = ", p)
  }

  nr <- n - p
  M <- matrix(0, nr, p)
  for (i in seq_len(p)) {
    M[, i] <- x[(p - i + 1):(n - i)]
  }
  M
}

#' Get Yule-Walker estimates for AR(p)
#'
#' @param x Numeric vector
#' @param p Integer value of AR order
#' @return Numeric vector of length p (AR coefficients)
#' @keywords internal
get_yw_estimates <- function(x, p) {
  # This is a simplified approach that may not handle edge cases
  r <- numeric(p+1)
  n <- length(x)
  xm <- mean(x)
  for (k in 0:p) {
    r[k+1] <- sum((x[1:(n-k)] - xm)*(x[(k+1):n] - xm))
  }
  R <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      R[i,j] <- r[abs(i-j)+1]
    }
  }
  rhs <- r[2:(p+1)]
  phi <- solve(R, rhs)
  phi
}

#' Create design matrix for time series
#'
#' @param x Time series data
#' @param model_info List with model parameters
#' @param innovations Optional innovations/residuals for MA components
#'
#' @return List with design matrix, response variable, and other components
#' @keywords internal
create_ts_design_matrix <- function(x, model_info, innovations = NULL) {
  # Extract model parameters
  ar_order <- model_info$ar_order
  ma_order <- model_info$ma_order
  d <- model_info$d
  model_type <- model_info$model_type
  include_mean <- model_info$include.mean

  # Check if series needs to be differenced
  if (model_type == "arima" && d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Handle mean
  if (include_mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  # Calculate maximum lag and effective data length
  max_lag <- max(ar_order, ma_order)
  n_data <- length(x_centered)

  if (n_data <= max_lag) {
    stop("Insufficient data to build model after differencing")
  }

  n_rows <- n_data - max_lag
  n_cols <- ar_order + ma_order

  # Create design matrix with appropriate dimensions
  X <- matrix(0, nrow = n_rows, ncol = n_cols)
  y <- x_centered[(max_lag + 1):n_data]

  # Add AR components
  if (ar_order > 0) {
    col_index <- 1
    for (i in 1:ar_order) {
      X[, col_index] <- x_centered[(max_lag - i + 1):(n_data - i)]
      col_index <- col_index + 1
    }
  }

  # Add MA components if needed
  if (ma_order > 0) {
    # If innovations not provided, use zeros
    if (is.null(innovations)) {
      innovations <- rep(0, n_data)
    }

    # Ensure correct length of innovations
    if (length(innovations) < n_data) {
      innovations <- c(rep(0, n_data - length(innovations)), innovations)
    } else if (length(innovations) > n_data) {
      innovations <- tail(innovations, n_data)
    }

    # Fill MA columns
    col_index <- ar_order + 1
    for (j in 1:ma_order) {
      X[, col_index] <- innovations[(max_lag - j + 1):(n_data - j)]
      col_index <- col_index + 1
    }
  }

  # Return results as a list
  list(
    X = X,
    y = y,
    x_centered = x_centered,
    x_mean = x_mean,
    n_rows = n_rows,
    n_cols = n_cols,
    effective_length = n_rows,
    innovations = innovations,
    original_x = x
  )
}

#' Get initial parameter estimates for time series models
#'
#' @param model_params Validated model parameters from validate_ts_parameters
#' @param initial Optionally user-provided initial estimates
#' @param method Estimation method
#' @param verbose Output detailed information
#'
#' @return List containing:
#'   \item{b_init}{Vector of initial AR/MA coefficients}
#'   \item{x_mean}{Estimated mean (if include.mean=TRUE)}
#'   \item{innovations}{Initial residuals/innovations}
#'   \item{x_centered}{Centered (or differenced + centered) series}
#'   \item{m2}{Second central moment of initial residuals}
#'   \item{m3}{Third central moment of initial residuals}
#'   \item{m4}{Fourth central moment of initial residuals}
#' @keywords internal
get_initial_estimates <- function(model_params,
                                  initial = NULL,
                                  method = "pmm2",
                                  verbose = FALSE) {
  x         <- model_params$original_x
  ar_order  <- model_params$ar_order
  ma_order  <- model_params$ma_order
  d         <- model_params$d
  mtype     <- model_params$model_type
  inc_mean  <- model_params$include.mean

  # Possibly difference for ARIMA
  if (mtype == "arima" && d > 0) {
    x_diff <- diff(x, differences = d)
  } else {
    x_diff <- x
  }

  # Center if needed
  if (inc_mean) {
    x_mean <- mean(x_diff, na.rm = TRUE)
    x_centered <- x_diff - x_mean
  } else {
    x_mean <- 0
    x_centered <- x_diff
  }

  if (mtype == "ar") {
    # AR(p): quick approach for initial values
    if (is.null(initial)) {
      if (method == "yw") {
        b_init <- get_yw_estimates(x_centered, ar_order)
      } else {
        X <- create_ar_matrix(x_centered, ar_order)
        y <- x_centered[(ar_order + 1):length(x_centered)]
        fit_ols <- lm.fit(x = X, y = y)
        b_init <- fit_ols$coefficients
      }
    } else {
      if (length(initial) != ar_order) {
        stop("Length of 'initial' must match AR order")
      }
      b_init <- initial
    }
    # Innovations from initial fit
    X <- create_ar_matrix(x_centered, ar_order)
    y <- x_centered[(ar_order + 1):length(x_centered)]
    innovations <- as.numeric(y - X %*% b_init)

  } else if (mtype %in% c("ma", "arma", "arima")) {
    # Use stats::arima for initial guess or user-provided
    arima_order <- c(ar_order, ifelse(mtype=="arima", d, 0), ma_order)
    if (is.null(initial)) {
      init_fit <- NULL
      try_methods <- c("CSS-ML","ML","CSS")

      for(mm in try_methods) {
        tmp <- tryCatch({
          stats::arima(x, order=arima_order, method=mm,
                       include.mean=inc_mean && (mtype!="arima"))
        }, error=function(e) NULL)
        if(!is.null(tmp)) {
          init_fit <- tmp
          break
        }
      }
      if(is.null(init_fit)) {
        if(verbose) cat("All standard methods failed; using simplified values.\n")
        init_fit <- list(
          coef = numeric(ar_order + ma_order),
          residuals = if(mtype=="arima") x_diff else x_centered
        )
        if(ar_order>0) names(init_fit$coef)[1:ar_order] <- paste0("ar",1:ar_order)
        if(ma_order>0) names(init_fit$coef)[(ar_order+1):(ar_order+ma_order)] <- paste0("ma",1:ma_order)
      }

      ar_init <- rep(0, ar_order)
      ma_init <- rep(0, ma_order)
      if(ar_order>0) {
        idx <- paste0("ar",1:ar_order)
        ar_init <- if(all(idx %in% names(init_fit$coef))) as.numeric(init_fit$coef[idx]) else rep(0.1, ar_order)
      }
      if(ma_order>0) {
        idx <- paste0("ma",1:ma_order)
        ma_init <- if(all(idx %in% names(init_fit$coef))) as.numeric(init_fit$coef[idx]) else rep(0.1, ma_order)
      }
      if(inc_mean && !is.null(init_fit$coef) && ("intercept" %in% names(init_fit$coef))) {
        x_mean <- init_fit$coef["intercept"]
      }
      innovations <- if(!is.null(init_fit$residuals)) as.numeric(init_fit$residuals) else {
        rnorm(length(x_centered),0,sd(x_centered,na.rm=TRUE))
      }
      b_init <- c(ar_init, ma_init)

    } else {
      # initial provided
      if(is.list(initial)) {
        if(ar_order>0 && is.null(initial$ar)) {
          stop("Missing 'ar' in initial list, but ar_order>0")
        }
        if(ma_order>0 && is.null(initial$ma)) {
          stop("Missing 'ma' in initial list, but ma_order>0")
        }
        ar_init <- if(ar_order>0) initial$ar else numeric(0)
        ma_init <- if(ma_order>0) initial$ma else numeric(0)
      } else {
        if(length(initial) != (ar_order+ma_order)) {
          stop("Length of 'initial' must match sum of AR and MA orders")
        }
        ar_init <- if(ar_order>0) initial[1:ar_order] else numeric(0)
        ma_init <- if(ma_order>0) initial[(ar_order+1):(ar_order+ma_order)] else numeric(0)
      }
      b_init <- c(ar_init, ma_init)

      init_fit <- tryCatch({
        stats::arima(x, order=arima_order,
                     fixed=b_init,
                     include.mean=inc_mean && (mtype!="arima"))
      }, error=function(e) {
        if(verbose) cat("Error with user-provided initial values:",e$message,"\n")
        list(residuals = if(mtype=="arima") x_diff else x_centered)
      })
      innovations <- as.numeric(init_fit$residuals)
    }
  }

  if(anyNA(b_init)) {
    warning("NA in initial parameters replaced with 0.")
    b_init[is.na(b_init)] <- 0
  }

  # Compute moments
  moments <- compute_moments(innovations)

  # Return results
  list(
    b_init      = b_init,
    x_mean      = x_mean,
    innovations = innovations,
    x_centered  = x_centered,
    orig_x      = x,
    m2          = moments$m2,
    m3          = moments$m3,
    m4          = moments$m4
  )
}

#' Update MA model innovations
#'
#' @param x Centered time series
#' @param ma_coef Vector of MA coefficients
#' @return Vector of innovations
#' @keywords internal
update_ma_innovations <- function(x, ma_coef) {
  n <- length(x)
  q <- length(ma_coef)

  # Initialize innovations as zeros
  innovations <- numeric(n)

  # Iteratively compute innovations
  for(t in 1:n) {
    # Calculate expected value based on previous innovations
    expected <- 0
    for(j in 1:q) {
      if(t - j > 0) {
        expected <- expected + ma_coef[j] * innovations[t - j]
      }
    }

    # Calculate current innovation
    innovations[t] <- x[t] - expected
  }

  # Check for infinite values
  if(any(is.infinite(innovations)) || any(is.na(innovations))) {
    warning("Infinite innovations detected in update_ma_innovations. Using regularized values.")
    # Replace problematic values with mean or 0
    bad_idx <- is.infinite(innovations) | is.na(innovations)
    if(sum(!bad_idx) > 0) {
      # If there are valid values, use their mean
      innovations[bad_idx] <- mean(innovations[!bad_idx])
    } else {
      # Otherwise use 0
      innovations[bad_idx] <- 0
    }
  }

  return(innovations)
}

#' Compute final residuals for time series models
#'
#' @param coefs Estimated coefficients
#' @param model_info Model information
#' @return Vector of residuals
#' @keywords internal
compute_ts_residuals <- function(coefs, model_info) {
  # Extract model parameters
  x           <- model_info$x
  ar_order    <- model_info$ar_order
  ma_order    <- model_info$ma_order
  d           <- model_info$d
  model_type  <- model_info$model_type
  include.mean <- model_info$include.mean
  verbose_flag <- isTRUE(model_info$verbose)

  # Separate coefficients into AR and MA parts
  if (ar_order > 0) {
    ar_coefs <- coefs[1:ar_order]
  } else {
    ar_coefs <- numeric(0)
  }

  if (ma_order > 0) {
    ma_coefs <- coefs[(ar_order+1):(ar_order+ma_order)]
  } else {
    ma_coefs <- numeric(0)
  }

  # For more reliable residual computation, use arima with fixed parameters
  arima_order <- c(ar_order, ifelse(model_type == "arima", d, 0), ma_order)

  # Prepare fixed parameters for arima
  fixed_params <- c(ar_coefs, ma_coefs)
  include_intercept <- include.mean && !(model_info$model_type == "arima" && model_info$d > 0)

  if (include_intercept) {
    fixed_params <- c(fixed_params, model_info$x_mean)
    names(fixed_params) <- c(
      if(ar_order > 0) paste0("ar", 1:ar_order) else NULL,
      if(ma_order > 0) paste0("ma", 1:ma_order) else NULL,
      "intercept"
    )
  } else {
    names(fixed_params) <- c(
      if(ar_order > 0) paste0("ar", 1:ar_order) else NULL,
      if(ma_order > 0) paste0("ma", 1:ma_order) else NULL
    )
  }

  # Compute residuals with fixed parameters
  final_fit <- tryCatch({
    stats::arima(model_info$original_x,
                order = arima_order,
                fixed = fixed_params,
                include.mean = include_intercept)
  }, error = function(e) {
    if (verbose_flag) cat("Error computing final residuals:", e$message, "\n")
    list(residuals = rep(NA, length(model_info$original_x)))
  })

  # Extract residuals and ensure correct length
  final_res <- as.numeric(final_fit$residuals)

  if (length(final_res) < length(model_info$original_x)) {
    final_res <- c(rep(NA, length(model_info$original_x) - length(final_res)), final_res)
  }

  return(final_res)
}
