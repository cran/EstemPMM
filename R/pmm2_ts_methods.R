# pmm2_ts_methods.R - Methods for working with time series model objects

#' Extract coefficients from TS2fit object
#'
#' @param object TS2fit object
#' @param ... Additional arguments (not used)
#'
#' @return Named vector of coefficients
#' @export
setMethod("coef", "TS2fit",
          function(object, ...) {
            # Get model parameters
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma

            # Extract and name AR coefficients
            if(ar_order > 0) {
              ar_coefs <- object@coefficients[1:ar_order]
              names(ar_coefs) <- paste0("ar", 1:ar_order)
            } else {
              ar_coefs <- numeric(0)
            }

            # Extract and name MA coefficients
            if(ma_order > 0) {
              ma_coefs <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
              names(ma_coefs) <- paste0("ma", 1:ma_order)
            } else {
              ma_coefs <- numeric(0)
            }

            # Combine coefficients
            result <- c(ar_coefs, ma_coefs)

            # Add intercept if present
            if(object@intercept != 0) {
              result <- c(intercept = object@intercept, result)
            }

            return(result)
          })

#' Extract residuals from TS2fit object
#'
#' @param object TS2fit object
#' @param ... Additional arguments (not used)
#'
#' @return Vector of residuals (innovations)
#' @export
setMethod("residuals", "TS2fit",
          function(object, ...) {
            object@residuals
          })

#' Get fitted values for AR model
#'
#' @param object TS2fit object with model_type="ar"
#' @return Vector of fitted values
#' @keywords internal
get_ar_fitted <- function(object) {
  if(object@model_type != "ar") {
    stop("This function is only for AR models")
  }

  x <- object@original_series
  ar_order <- object@order$ar
  ar_coef <- object@coefficients[1:ar_order]
  intercept <- object@intercept

  if(intercept != 0) {
    x_centered <- x - intercept
  } else {
    x_centered <- x
  }

  # Create design matrix and calculate fitted values
  X <- create_ar_matrix(x_centered, ar_order)
  fitted <- as.vector(X %*% ar_coef) + intercept
  return(fitted)
}

#' Extract fitted values from TS2fit object
#'
#' @param object TS2fit object
#' @param ... Additional arguments (not used)
#'
#' @return Vector of fitted values
#' @export
setMethod("fitted", "TS2fit",
          function(object, ...) {
            # Get model type
            model_type <- object@model_type

            # Calculate fitted values based on model type
            if(model_type == "ar") {
              # For AR models use direct calculation
              fitted_values <- get_ar_fitted(object)
            } else {
              # For other models: fitted = original minus residuals
              orig <- object@original_series
              resid <- object@residuals

              # Align lengths (often residuals are shorter due to initial values)
              len_diff <- length(orig) - length(resid)
              if(len_diff > 0 && !all(is.na(resid))) {
                # Find first non-NA value in resid
                first_valid <- min(which(!is.na(resid)))

                # Build fitted vector with NA at the beginning
                fitted_values <- rep(NA, length(orig))
                valid_indices <- first_valid:length(resid)

                # Set actual values
                fitted_values[(len_diff + valid_indices)] <-
                  orig[(len_diff + valid_indices)] - resid[valid_indices]
              } else {
                fitted_values <- orig - resid
              }
            }

            return(fitted_values)
          })

#' Build diagnostic plots for TS2fit objects
#'
#' @param x TS2fit object
#' @param y Not used (for S4 method compatibility)
#' @param which Integer vector indicating which plots to produce
#' @param ... additional arguments passed to plot functions
#'
#' @return Invisibly returns x
#'
#' @export
setMethod("plot", signature(x = "TS2fit", y = "missing"),
          function(x, y, which = c(1:4), ...) {
            op <- par(no.readonly = TRUE)
            on.exit(par(op))

            # Get model parameters
            model_type <- x@model_type
            ar_order <- x@order$ar
            ma_order <- x@order$ma
            d <- x@order$d

            # Default plot layout
            par(mfrow = c(2, 2))

            # For ARIMA models we may want to plot original/differenced series as well
            if(model_type == "arima" && length(which) > 4) {
              par(mfrow = c(3, 2))
            }

            # Calculate fitted values and residuals
            residuals <- as.numeric(x@residuals)
            fitted <- fitted(x)

            # Determine which plots to display
            plot_idx <- 1
            n_plots <- min(length(which), 6) # Maximum 6 plots

            # For ARIMA models, we may want other plots
            if(model_type == "arima") {
              # Plot 1: Original time series (ARIMA only)
              if(1 %in% which && plot_idx <= n_plots) {
                plot(x@original_series, type = "l",
                     main = "Original Time Series",
                     xlab = "Time",
                     ylab = "Value",
                     ...)
                plot_idx <- plot_idx + 1
              }

              # Plot 2: Differenced time series (ARIMA only)
              if(2 %in% which && plot_idx <= n_plots && d > 0) {
                diff_series <- diff(x@original_series, differences = d)
                plot(diff_series, type = "l",
                     main = paste0("Differenced Series (d=", d, ")"),
                     xlab = "Time",
                     ylab = "Value",
                     ...)
                plot_idx <- plot_idx + 1
              }
            }

            # Standard plots for all model types
            # Plot: Residuals vs Fitted
            if(3 %in% which && plot_idx <= n_plots) {
              plot(fitted, residuals,
                   main = "Residuals vs Fitted",
                   xlab = "Fitted values",
                   ylab = "Residuals",
                   ...)
              abline(h = 0, lty = 2)
              lines(lowess(fitted, residuals), col = "red")
              plot_idx <- plot_idx + 1
            }

            # Plot: Normal Q-Q plot
            if(4 %in% which && plot_idx <= n_plots) {
              qqnorm(residuals, main = "Normal Q-Q Plot", ...)
              qqline(residuals)
              plot_idx <- plot_idx + 1
            }

            # Plot: ACF of residuals
            if(5 %in% which && plot_idx <= n_plots) {
              acf(residuals, main = "ACF of Residuals", ...)
              plot_idx <- plot_idx + 1
            }

            # Plot: Histogram of residuals
            if(6 %in% which && plot_idx <= n_plots) {
              hist(residuals,
                   main = "Histogram of Residuals",
                   xlab = "Residuals",
                   breaks = "FD",
                   ...)
              plot_idx <- plot_idx + 1
            }

            invisible(x)
          })

#' Prediction method for TS2fit objects
#'
#' @param object TS2fit object
#' @param n.ahead Number of steps ahead for prediction
#' @param ... additional arguments (not used)
#'
#' @return Vector or list of predictions, depending on model type
#'
#' @export
setMethod("predict", "TS2fit",
          function(object, n.ahead = 1, ...) {
            # Get model parameters
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d
            intercept <- object@intercept

            # Extract coefficients
            if(ar_order > 0) {
              ar_coef <- object@coefficients[1:ar_order]
            } else {
              ar_coef <- numeric(0)
            }

            if(ma_order > 0) {
              ma_coef <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
            } else {
              ma_coef <- numeric(0)
            }

            # For AR models, implement direct prediction
            if(model_type == "ar") {
              x <- object@original_series
              n <- length(x)
              pred <- numeric(n.ahead)

              # Generate predictions
              for(i in 1:n.ahead) {
                # Use original data and previous predictions as needed
                lags <- numeric(ar_order)
                for(j in 1:ar_order) {
                  if(i - j <= 0) {
                    # Use original data
                    lags[j] <- x[n - j + i]
                  } else {
                    # Use previous predictions
                    lags[j] <- pred[i - j]
                  }
                }

                # Calculate prediction
                pred[i] <- sum(ar_coef * lags) + intercept
              }

              return(pred)

            } else if(model_type == "ma") {
              # For MA models, predictions beyond the order are just the mean
              innovations <- object@residuals

              # Generate MA predictions
              ma_pred <- function(innovations, ma_coef, n.ahead) {
                n <- length(innovations)
                q <- length(ma_coef)
                pred <- numeric(n.ahead)

                for(i in 1:n.ahead) {
                  for(j in 1:min(i, q)) {
                    if((n - i + j) > 0) {
                      pred[i] <- pred[i] + ma_coef[j] * innovations[n - i + j]
                    }
                  }
                }
                return(pred)
              }

              if(n.ahead > ma_order) {
                return(c(ma_pred(innovations, ma_coef, ma_order),
                         rep(intercept, n.ahead - ma_order)))
              } else {
                return(ma_pred(innovations, ma_coef, n.ahead))
              }

            } else {
              # For ARMA and ARIMA models, use stats::arima predictions
              # which correctly handle both components

              # Set up arima model with fixed parameters
              arima_order <- c(ar_order, ifelse(model_type == "arima", d, 0), ma_order)

              # Use prediction function from stats package
              arima_pred <- stats::predict(
                stats::arima(object@original_series,
                             order = arima_order,
                             include.mean = (intercept != 0),
                             fixed = c(ar_coef, ma_coef, if(intercept != 0) intercept else NULL)),
                n.ahead = n.ahead
              )

              return(arima_pred)
            }
          })

#' Compare PMM2 with classical time series estimation methods
#'
#' @param x Numeric vector of time series data
#' @param order Model order specification (see ts_pmm2 for format)
#' @param model_type Model type: "ar", "ma", "arma", or "arima"
#' @param include.mean Logical, whether to include intercept term
#' @param pmm2_args List of additional arguments to pass to ts_pmm2()
#'
#' @return A named list containing the fitted objects for each estimation
#'   approach (e.g., YW/OLS/MLE or CSS/ML alongside PMM2) plus two data frames:
#'   \code{coefficients} (side-by-side parameter estimates) and
#'   \code{residual_stats} (residual RSS, MAE, skewness, and kurtosis).
#' @export
compare_ts_methods <- function(x, order, model_type = c("ar", "ma", "arma", "arima"),
                               include.mean = TRUE, pmm2_args = list()) {
  # Select model_type argument
  model_type <- match.arg(model_type)

  # Prepare model comparison based on model_type
  if(model_type == "ar") {
    # For AR models
    # Fit AR model using Yule-Walker method
    yw_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "yw",
                        demean = include.mean)

    # Fit AR model using OLS method
    ols_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "ols",
                         demean = include.mean)

    # Fit AR model using MLE method
    mle_fit <- stats::ar(x, order.max = order, aic = FALSE, method = "mle",
                         demean = include.mean)

    # Fit AR model using PMM2
    pmm2_args <- c(list(x = x, order = order, model_type = "ar",
                        include.mean = include.mean), pmm2_args)
    pmm2_fit <- do.call(ts_pmm2, pmm2_args)

    # Extract coefficients
    coef_yw <- yw_fit$ar
    coef_ols <- ols_fit$ar
    coef_mle <- mle_fit$ar
    coef_pmm2 <- pmm2_fit@coefficients

    # Calculate residuals
    res_yw <- yw_fit$resid[!is.na(yw_fit$resid)]
    res_ols <- ols_fit$resid[!is.na(ols_fit$resid)]
    res_mle <- mle_fit$resid[!is.na(mle_fit$resid)]
    res_pmm2 <- pmm2_fit@residuals

    methods <- c("YW", "OLS", "MLE", "PMM2")

    result_list <- list(
      yw = yw_fit,
      ols = ols_fit,
      mle = mle_fit,
      pmm2 = pmm2_fit
    )

  } else if(model_type %in% c("ma", "arma", "arima")) {
    # For MA, ARMA and ARIMA models

    # Prepare arima order based on model type
    if(model_type == "ma") {
      arima_order <- c(0, 0, order)
    } else if(model_type == "arma") {
      arima_order <- c(order[1], 0, order[2])
    } else {
      arima_order <- order
    }

    # Fit model using CSS method
    css_fit <- arima(x, order = arima_order, method = "CSS", include.mean = include.mean)

    # Fit model using ML method
    ml_fit <- arima(x, order = arima_order, method = "ML", include.mean = include.mean)

    # Fit model using PMM2
    pmm2_args <- c(list(x = x, order = order, model_type = model_type,
                        include.mean = include.mean), pmm2_args)
    pmm2_fit <- do.call(ts_pmm2, pmm2_args)

    # Extract AR and MA coefficient names based on model type
    if(model_type == "ma") {
      ar_names <- character(0)
      ma_names <- paste0("ma", 1:order)
    } else if(model_type == "arma") {
      ar_names <- paste0("ar", 1:order[1])
      ma_names <- paste0("ma", 1:order[2])
    } else {
      ar_names <- if(order[1] > 0) paste0("ar", 1:order[1]) else character(0)
      ma_names <- if(order[3] > 0) paste0("ma", 1:order[3]) else character(0)
    }

    coef_names <- c(ar_names, ma_names)

    # Extract coefficients
    coef_css <- as.numeric(css_fit$coef[coef_names])
    coef_ml <- as.numeric(ml_fit$coef[coef_names])
    coef_pmm2 <- pmm2_fit@coefficients

    # Calculate residuals
    res_css <- residuals(css_fit)
    res_ml <- residuals(ml_fit)
    res_pmm2 <- pmm2_fit@residuals

    methods <- c("CSS", "ML", "PMM2")

    result_list <- list(
      css = css_fit,
      ml = ml_fit,
      pmm2 = pmm2_fit
    )
  }

  # Calculate residual statistics for all methods
  residuals_list <- if(model_type == "ar") {
    list(res_yw, res_ols, res_mle, res_pmm2)
  } else {
    list(res_css, res_ml, res_pmm2)
  }

  compute_res_stats <- function(res) {
    m2 <- mean(res^2, na.rm = TRUE)
    m3 <- mean(res^3, na.rm = TRUE)
    m4 <- mean(res^4, na.rm = TRUE)

    c(RSS = sum(res^2, na.rm = TRUE),
      MAE = mean(abs(res), na.rm = TRUE),
      Skewness = m3 / m2^(3/2),
      Kurtosis = m4 / m2^2)
  }

  res_stats <- data.frame(
    Method = methods,
    do.call(rbind, lapply(residuals_list, compute_res_stats))
  )

  # Create coefficient comparison table
  if(model_type == "ar") {
    coef_names <- paste0("ar", 1:order)
    coef_values <- list(coef_yw, coef_ols, coef_mle, coef_pmm2)
  } else if(model_type == "ma") {
    coef_names <- paste0("ma", 1:order)
    coef_values <- list(coef_css, coef_ml, coef_pmm2)
  } else {
    coef_values <- list(coef_css, coef_ml, coef_pmm2)
  }

  coef_table <- data.frame(
    Coefficient = coef_names,
    do.call(cbind, lapply(seq_along(methods), function(i) {
      result <- coef_values[[i]]
      names(result) <- methods[i]
      return(result)
    }))
  )

  # Return results
  result_list$coefficients <- coef_table
  result_list$residual_stats <- res_stats

  return(result_list)
}

#' Compare AR methods
#'
#' @inheritParams compare_ts_methods
#' @inherit compare_ts_methods return
#' @export
compare_ar_methods <- function(x, order = 1, include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "ar",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Compare MA methods
#'
#' @inheritParams compare_ts_methods
#' @inherit compare_ts_methods return
#' @export
compare_ma_methods <- function(x, order = 1, include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "ma",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Compare ARMA methods
#'
#' @inheritParams compare_ts_methods
#' @inherit compare_ts_methods return
#' @export
compare_arma_methods <- function(x, order = c(1, 1), include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "arma",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}

#' Compare ARIMA methods
#'
#' @inheritParams compare_ts_methods
#' @inherit compare_ts_methods return
#' @export
compare_arima_methods <- function(x, order = c(1, 1, 1), include.mean = TRUE, pmm2_args = list()) {
  compare_ts_methods(x, order = order, model_type = "arima",
                     include.mean = include.mean, pmm2_args = pmm2_args)
}
