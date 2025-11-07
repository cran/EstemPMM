# pmm2_inference.R - Statistical inference for PMM2 models

#' Bootstrap inference for PMM2 fit
#'
#' @param object object of class PMM2fit
#' @param formula the same formula that was used initially
#' @param data data frame that was used initially
#' @param B number of bootstrap replications
#' @param seed (optional) for reproducibility
#' @param parallel logical, whether to use parallel computing
#' @param cores number of cores to use for parallel computing, defaults to auto-detect
#'
#' @return data.frame with columns: Estimate, Std.Error, t.value, p.value
#' @export
pmm2_inference <- function(object, formula, data, B=200, seed=NULL,
                           parallel=FALSE, cores=NULL) {
  # Set seed for reproducibility if provided
  if(!is.null(seed)) set.seed(seed)

  # Extract coefficients and residuals
  coefs <- object@coefficients
  res   <- object@residuals

  # Validate input data
  if(B < 10) {
    warning("Number of bootstrap samples (B) is too small. Consider using B >= 100 for more reliable inference.")
  }

  if(!inherits(object, "PMM2fit")) {
    stop("Object must be of class 'PMM2fit'")
  }

  if(missing(formula) || missing(data)) {
    stop("Both 'formula' and 'data' must be provided")
  }

  # Build matrices X, y
  mf <- model.frame(formula, data)
  X <- model.matrix(formula, mf)
  y <- model.response(mf)
  n <- nrow(X)

  # Early return in case of errors
  if(is.null(y) || is.null(X)) {
    stop("Failed to extract response or design matrix from data")
  }

  # Check whether to use parallel computing
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE)

  if(use_parallel) {
    if(is.null(cores)) {
      cores <- max(1, parallel::detectCores() - 1)
    }

    boot_results <- parallel::mclapply(seq_len(B), function(b) {
      # 1) Bootstrap residuals
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Create new y
      y_b <- X %*% coefs + res_b

      # 3) Create new data
      data_b <- data
      # Assume left side is the first term in formula
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Re-estimate model
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Bootstrap replication ", b, " failed: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }, mc.cores = cores)

    # Convert list to matrix
    boot_est <- do.call(rbind, boot_results)

  } else {
    # Sequential computing
    # Matrix to store results
    boot_est <- matrix(0, nrow=B, ncol=length(coefs))
    colnames(boot_est) <- names(coefs)

    # Progress tracking
    pb <- NULL
    if(interactive() && B > 10) {
      if(requireNamespace("utils", quietly = TRUE)) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    for(b in seq_len(B)) {
      # 1) Bootstrap residuals
      res_b <- sample(res, size=n, replace=TRUE)

      # 2) Create new y
      y_b <- X %*% coefs + res_b

      # 3) Create new data
      data_b <- data
      # Assume left side is the first term in formula
      lhs <- as.character(formula[[2]])
      data_b[[lhs]] <- as.numeric(y_b)

      # 4) Re-estimate model
      fit_b <- tryCatch({
        lm_pmm2(formula, data_b, max_iter=20, tol=1e-6)
      }, error = function(e) {
        warning("Bootstrap replication ", b, " failed: ", e$message)
        return(NULL)
      })

      if(!is.null(fit_b)) {
        boot_est[b, ] <- fit_b@coefficients
      } else {
        boot_est[b, ] <- NA
      }

      # Update progress indicator
      if(!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }

    # Close progress indicator
    if(!is.null(pb)) close(pb)
  }

  # Remove rows with NA values
  na_rows <- apply(boot_est, 1, function(row) any(is.na(row)))
  if(any(na_rows)) {
    warning("Removed ", sum(na_rows), " bootstrap replications due to estimation errors")
    boot_est <- boot_est[!na_rows, , drop = FALSE]
  }

  # Check that we have enough successful bootstraps
  if(nrow(boot_est) < 10) {
    stop("Too few successful bootstrap replications to compute reliable inference")
  }

  # Compute covariance matrix and standard errors
  cov_mat <- cov(boot_est)
  est <- coefs
  se  <- sqrt(diag(cov_mat))

  # Compute t-values and p-values
  t_val <- est / se
  # For large samples use normal approximation
  p_val <- 2 * (1 - pnorm(abs(t_val)))

  # Create output data frame
  out <- data.frame(
    Estimate  = est,
    Std.Error = se,
    t.value   = t_val,
    p.value   = p_val
  )
  rownames(out) <- names(est)

  # Compute confidence intervals
  ci <- t(apply(boot_est, 2, quantile, probs = c(0.025, 0.975)))
  colnames(ci) <- c("2.5%", "97.5%")

  # Add confidence intervals to output
  out$conf.low <- ci[, "2.5%"]
  out$conf.high <- ci[, "97.5%"]

  return(out)
}

#' Plot bootstrap distributions for PMM2 fit
#'
#' @param object Result from pmm2_inference
#' @param coefficients Which coefficients to plot, defaults to all
#'
#' @return Invisibly returns histogram information
#' @export
plot_pmm2_bootstrap <- function(object, coefficients = NULL) {
  if(!inherits(object, "data.frame") ||
     !all(c("Estimate", "Std.Error", "conf.low", "conf.high") %in% names(object))) {
    stop("Object must be the result from pmm2_inference()")
  }

  # If coefficients not specified, use all
  if(is.null(coefficients)) {
    coefficients <- rownames(object)
  }

  # Filter to requested coefficients
  object_subset <- object[intersect(coefficients, rownames(object)), , drop = FALSE]

  # Check for empty dataset
  if(nrow(object_subset) == 0) {
    warning("None of the requested coefficients found in results.")
    return(invisible(NULL))
  }

  # Set up plot layout
  n_coefs <- nrow(object_subset)
  n_cols <- min(2, n_coefs)
  n_rows <- ceiling(n_coefs / n_cols)

  # Save old par settings and restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  par(mfrow = c(n_rows, n_cols))

  # Create density plot for each coefficient
  result <- lapply(seq_len(n_coefs), function(i) {
    coef_name <- rownames(object_subset)[i]
    est <- object_subset[i, "Estimate"]
    ci_low <- object_subset[i, "conf.low"]
    ci_high <- object_subset[i, "conf.high"]
    se <- object_subset[i, "Std.Error"]

    # Check for finite values
    if(!is.finite(est) || !is.finite(ci_low) || !is.finite(ci_high) || !is.finite(se)) {
      warning("Infinite or NA values for coefficient ", coef_name,
              ". Skipping this plot.")
      return(NULL)
    }

    # Create plot title
    main_title <- paste0(coef_name, "\nEstimate: ", round(est, 4))

    # Estimate range for x-axis
    # Use more robust approach to determine range
    x_range <- range(c(est, ci_low, ci_high), na.rm = TRUE)
    # Expand range by 20% in both directions
    x_range_width <- diff(x_range)
    x_range <- x_range + c(-0.2, 0.2) * x_range_width

    # Create points for x-axis
    x_seq <- seq(x_range[1], x_range[2], length.out = 100)

    # Create density values for normal distribution
    y_seq <- dnorm(x_seq, mean = est, sd = se)

    # Plot
    plot(x_seq, y_seq, type = "l",
         main = main_title,
         xlab = "Value",
         ylab = "Density")

    # Add vertical lines for estimate and CI
    abline(v = est, col = "red", lwd = 2)
    abline(v = ci_low, col = "blue", lty = 2)
    abline(v = ci_high, col = "blue", lty = 2)

    # Add legend
    legend("topright",
           legend = c("Estimate", "95% CI"),
           col = c("red", "blue"),
           lty = c(1, 2),
           lwd = c(2, 1),
           cex = 0.8)

    invisible(list(x = x_seq, y = y_seq, estimate = est,
                   ci_low = ci_low, ci_high = ci_high))
  })

  # Remove NULL results
  result <- result[!sapply(result, is.null)]

  # If all results are NULL, return NULL
  if(length(result) == 0) {
    warning("Failed to create any plots.")
    return(invisible(NULL))
  }

  # Add names to results
  names(result) <- rownames(object_subset)[sapply(seq_len(n_coefs), function(i) {
    !is.null(result[[i]])
  })]

  invisible(result)
}


#' Bootstrap inference for PMM2 time series models
#'
#' @param object object of class TS2fit
#' @param x (optional) original time series; if NULL, uses object@original_series
#' @param B number of bootstrap replications
#' @param seed (optional) for reproducibility
#' @param block_length block length for block bootstrap; if NULL, uses heuristic value
#' @param method bootstrap type: "residual" or "block"
#' @param parallel logical, whether to use parallel computing
#' @param cores number of cores for parallel computing
#' @param debug logical, whether to output additional diagnostic information
#'
#' @return data.frame with columns: Estimate, Std.Error, t.value, p.value
#' @export
ts_pmm2_inference <- function(object, x = NULL, B = 200, seed = NULL,
                              block_length = NULL, method = c("residual", "block"),
                              parallel = FALSE, cores = NULL, debug = FALSE) {
  # Check object class
  if (!inherits(object, "TS2fit")) {
    stop("Object must be of class 'TS2fit'")
  }

  # Select bootstrap method
  method <- match.arg(method)

  # Set seed for reproducibility if provided
  if (!is.null(seed)) set.seed(seed)

  # Extract model parameters
  model_type <- object@model_type
  ar_order <- object@order$ar
  ma_order <- object@order$ma
  d <- object@order$d
  intercept <- object@intercept
  include_mean <- intercept != 0

  if(debug) {
    cat("Model parameters:\n")
    cat("model_type:", model_type, "\n")
    cat("ar_order:", ar_order, "\n")
    cat("ma_order:", ma_order, "\n")
    cat("d:", d, "\n")
    cat("intercept:", intercept, "\n")
    cat("include_mean:", include_mean, "\n")
  }

  # If x not provided, use original series from object
  if (is.null(x)) {
    x <- object@original_series
  }

  # Extract coefficients and residuals
  coefs <- object@coefficients
  res <- object@residuals

  if(debug) {
    cat("Original series size:", length(x), "\n")
    cat("Residuals vector size:", length(res), "\n")
    cat("Number of coefficients:", length(coefs), "\n")
  }

  # Perform block bootstrap if specified
  if (method == "block") {
    # Determine block length if not provided
    if (is.null(block_length)) {
      # Use heuristic: square root of series length
      block_length <- ceiling(sqrt(length(x)))
    }

    # Check that block length makes sense
    if (block_length < 2) {
      warning("Block length too small. Setting to 2.")
      block_length <- 2
    }
    if (block_length > length(x) / 4) {
      warning("Block length too large. Setting to 1/4 of series length.")
      block_length <- floor(length(x) / 4)
    }

    # Function to generate bootstrap series from block bootstrap
    generate_block_bootstrap <- function(x, block_length) {
      n <- length(x)
      blocks_needed <- ceiling(n / block_length)

      # Possible block start positions
      start_positions <- 1:(n - block_length + 1)

      # Select random start positions
      selected_starts <- sample(start_positions, blocks_needed, replace = TRUE)

      # Create bootstrap series
      boot_series <- numeric(0)
      for (start in selected_starts) {
        boot_series <- c(boot_series, x[start:(start + block_length - 1)])
      }

      # Trim to original length
      boot_series[1:n]
    }

    # Bootstrap logic differs for block method
    boot_function <- function(b) {
      # Generate new series using block bootstrap
      x_b <- generate_block_bootstrap(x, block_length)

      # Determine proper order format depending on model type
      if(model_type == "ar") {
        boot_order <- ar_order  # For AR models - single number
      } else if(model_type == "ma") {
        boot_order <- ma_order  # For MA models - single number
      } else if(model_type == "arma") {
        boot_order <- c(ar_order, ma_order)  # For ARMA - vector of length 2
      } else if(model_type == "arima") {
        boot_order <- c(ar_order, d, ma_order)  # For ARIMA - vector of length 3
      } else {
        stop("Unknown model type: ", model_type)
      }

      # Fit model on bootstrap series
      fit_b <- tryCatch({
        ts_pmm2(x_b, order = boot_order,
                model_type = model_type,
                include.mean = include_mean)
      }, error = function(e) {
        warning("Bootstrap replication ", b, " failed: ", e$message)
        return(NULL)
      })

      if (!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }
  } else {
    # For residual bootstrap
    # Function to generate new series based on model and bootstrap residuals
    generate_with_residuals <- function(model, residuals) {
      n <- length(model@original_series)

      # For ARIMA models, need to differentiate first
      if (model_type == "arima" && d > 0) {
        # Here need more complex logic to restore original series
        # after generating differenced series
        # This is a simplified approach:
        diff_x <- numeric(n - d)

        # Bootstrap residuals
        res_b <- sample(residuals[!is.na(residuals)], size = n - d, replace = TRUE)

        # Generate differenced series
        if (ar_order > 0) {
          ar_coefs <- model@coefficients[1:ar_order]
        } else {
          ar_coefs <- numeric(0)
        }

        if (ma_order > 0) {
          ma_coefs <- model@coefficients[(ar_order+1):(ar_order+ma_order)]
        } else {
          ma_coefs <- numeric(0)
        }

        # Simplified ARIMA process simulation
        diff_x <- arima.sim(model = list(
          ar = if(ar_order > 0) ar_coefs else NULL,
          ma = if(ma_order > 0) ma_coefs else NULL
        ), n = n - d, innov = res_b, n.start = max(ar_order, ma_order))

        # Integrate back
        x_b <- diffinv(diff_x, differences = d)

        # Add intercept if needed
        if (include_mean) {
          x_b <- x_b + intercept
        }
      } else {
        # For AR, MA, ARMA models
        res_b <- sample(residuals[!is.na(residuals)], size = n, replace = TRUE)

        # Extract AR and MA coefficients
        if (ar_order > 0) {
          ar_coefs <- model@coefficients[1:ar_order]
        } else {
          ar_coefs <- numeric(0)
        }

        if (ma_order > 0) {
          ma_coefs <- model@coefficients[(ar_order+1):(ar_order+ma_order)]
        } else {
          ma_coefs <- numeric(0)
        }

        # Simulate ARMA process
        x_b <- arima.sim(model = list(
          ar = if(ar_order > 0) ar_coefs else NULL,
          ma = if(ma_order > 0) ma_coefs else NULL
        ), n = n, innov = res_b, n.start = max(ar_order, ma_order))

        # Add intercept if needed
        if (include_mean) {
          x_b <- x_b + intercept
        }
      }

      return(x_b)
    }

    boot_function <- function(b) {
      # Generate new series
      x_b <- generate_with_residuals(object, res)

      # Determine proper order format depending on model type
      if(model_type == "ar") {
        boot_order <- ar_order  # For AR models - single number
      } else if(model_type == "ma") {
        boot_order <- ma_order  # For MA models - single number
      } else if(model_type == "arma") {
        boot_order <- c(ar_order, ma_order)  # For ARMA - vector of length 2
      } else if(model_type == "arima") {
        boot_order <- c(ar_order, d, ma_order)  # For ARIMA - vector of length 3
      } else {
        stop("Unknown model type: ", model_type)
      }

      if(debug && b == 1) {
        cat("Bootstrap replication 1:\n")
        cat("Model type:", model_type, "\n")
        cat("Order for bootstrap:", paste(boot_order, collapse=", "), "\n")
      }

      # Fit model on bootstrap series
      fit_b <- tryCatch({
        ts_pmm2(x_b, order = boot_order,
                model_type = model_type,
                include.mean = include_mean)
      }, error = function(e) {
        warning("Bootstrap replication ", b, " failed: ", e$message)
        return(NULL)
      })

      if (!is.null(fit_b)) {
        return(fit_b@coefficients)
      } else {
        return(rep(NA, length(coefs)))
      }
    }
  }

  # Perform bootstrap: parallel or sequential
  use_parallel <- parallel && requireNamespace("parallel", quietly = TRUE)

  if (use_parallel) {
    if (is.null(cores)) {
      cores <- max(1, parallel::detectCores() - 1)
    }

    boot_results <- parallel::mclapply(seq_len(B), function(b) {
      boot_function(b)
    }, mc.cores = cores)

    # Convert list to matrix
    boot_est <- do.call(rbind, boot_results)
  } else {
    # Sequential computing
    boot_est <- matrix(0, nrow = B, ncol = length(coefs))
    colnames(boot_est) <- names(coefs)

    # Progress tracking
    pb <- NULL
    if (interactive() && B > 10) {
      if (requireNamespace("utils", quietly = TRUE)) {
        pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
      }
    }

    for (b in seq_len(B)) {
      boot_est[b, ] <- boot_function(b)

      # Update progress indicator
      if (!is.null(pb)) utils::setTxtProgressBar(pb, b)
    }

    # Close progress indicator
    if (!is.null(pb)) close(pb)
  }

  # Remove rows with NA values
  na_rows <- apply(boot_est, 1, function(row) any(is.na(row)))
  if (any(na_rows)) {
    warning("Removed ", sum(na_rows), " bootstrap replications due to estimation errors")
    boot_est <- boot_est[!na_rows, , drop = FALSE]
  }

  # Check that we have enough successful bootstraps
  if (nrow(boot_est) < 10) {
    stop("Too few successful bootstrap replications to compute reliable inference")
  }

  # Check for NaN or Inf values
  if (any(is.nan(boot_est)) || any(is.infinite(boot_est))) {
    warning("Detected NaN or infinite values in bootstrap replications. Replacing them with NA.")
    boot_est[is.nan(boot_est) | is.infinite(boot_est)] <- NA
  }

  # Compute covariance matrix and standard errors
  cov_mat <- cov(boot_est, use = "pairwise.complete.obs")
  est <- coefs
  se <- sqrt(diag(cov_mat))

  # Compute t-values and p-values
  t_val <- est / se
  p_val <- 2 * (1 - pnorm(abs(t_val)))

  # Create output data frame
  out <- data.frame(
    Estimate = est,
    Std.Error = se,
    t.value = t_val,
    p.value = p_val
  )

  # Add names for AR and MA parameters
  param_names <- c()
  if (ar_order > 0) {
    param_names <- c(param_names, paste0("ar", 1:ar_order))
  }
  if (ma_order > 0) {
    param_names <- c(param_names, paste0("ma", 1:ma_order))
  }
  rownames(out) <- param_names

  # Compute confidence intervals
  ci <- t(apply(boot_est, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  colnames(ci) <- c("2.5%", "97.5%")

  # Add confidence intervals to output
  out$conf.low <- ci[, "2.5%"]
  out$conf.high <- ci[, "97.5%"]

  return(out)
}
