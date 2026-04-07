# pmm3_ts_methods.R - S4 methods for PMM3 time series objects

#' Extract coefficients from TS3fit object
#'
#' @param object TS3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Named vector of coefficients
#' @export
setMethod("coef", "TS3fit",
          function(object, ...) {
            ar_order <- object@order$ar
            ma_order <- object@order$ma

            ar_coefs <- if (ar_order > 0) {
              v <- object@coefficients[1:ar_order]
              names(v) <- paste0("ar", 1:ar_order)
              v
            } else {
              numeric(0)
            }

            ma_coefs <- if (ma_order > 0) {
              v <- object@coefficients[(ar_order + 1):(ar_order + ma_order)]
              names(v) <- paste0("ma", 1:ma_order)
              v
            } else {
              numeric(0)
            }

            result <- c(ar_coefs, ma_coefs)
            if (object@intercept != 0) {
              result <- c(intercept = object@intercept, result)
            }
            result
          })

#' Extract residuals from TS3fit object
#'
#' @param object TS3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Vector of residuals
#' @export
setMethod("residuals", "TS3fit",
          function(object, ...) {
            object@residuals
          })

#' Extract fitted values from TS3fit object
#'
#' @param object TS3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Vector of fitted values
#' @export
setMethod("fitted", "TS3fit",
          function(object, ...) {
            x <- object@original_series
            res <- object@residuals

            if (object@model_type == "ar") {
              ar_order <- object@order$ar
              ar_coef  <- object@coefficients[1:ar_order]
              intercept <- object@intercept

              x_centered <- if (intercept != 0) x - intercept else x
              X <- create_ar_matrix(x_centered, ar_order)
              fitted_vals <- as.vector(X %*% ar_coef) + intercept
              return(fitted_vals)
            }

            # Fallback: original - residuals (align lengths)
            n_res <- length(res)
            n_x   <- length(x)
            if (n_res < n_x) {
              return(x[(n_x - n_res + 1):n_x] - res)
            }
            x - res
          })

#' Summary method for TS3fit objects
#'
#' @param object TS3fit object
#' @param ... Additional arguments (not used)
#'
#' @return Prints summary to console; returns object (invisibly)
#' @export
setMethod("summary", "TS3fit",
          function(object, ...) {
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d

            cat("PMM3 time series estimation results (S = 3, symmetric errors)\n")

            cat("Model type: ")
            if (model_type == "ar") {
              cat("AR(", ar_order, ")\n", sep = "")
            } else if (model_type == "ma") {
              cat("MA(", ma_order, ")\n", sep = "")
            } else if (model_type == "arma") {
              cat("ARMA(", ar_order, ",", ma_order, ")\n", sep = "")
            } else if (model_type == "arima") {
              cat("ARIMA(", ar_order, ",", d, ",", ma_order, ")\n", sep = "")
            }

            if (!is.null(object@call)) {
              cat("Call:\n")
              print(object@call)
              cat("\n")
            }

            cat("Coefficients:\n")
            if (ar_order > 0) {
              ar_coefs <- object@coefficients[1:ar_order]
              names(ar_coefs) <- paste0("ar", 1:ar_order)
              cat("AR: ")
              print(ar_coefs)
            }
            if (ma_order > 0) {
              ma_coefs <- object@coefficients[(ar_order + 1):(ar_order + ma_order)]
              names(ma_coefs) <- paste0("ma", 1:ma_order)
              cat("MA: ")
              print(ma_coefs)
            }
            if (object@intercept != 0) {
              cat("Intercept: ", object@intercept, "\n")
            }

            cat("\nCentral moments of initial residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m4 =", object@m4, "\n")
            cat("  m6 =", object@m6, "\n\n")

            cat("Theoretical characteristics of PMM3 (S = 3):\n")
            cat("  gamma4 =", object@gamma4, "\n")
            cat("  gamma6 =", object@gamma6, "\n")
            cat("  g3     =", object@g_coefficient,
                " (expected ratio Var[PMM3]/Var[OLS])\n")
            if (!is.na(object@kappa)) {
              cat("  kappa  =", object@kappa, "\n")
            }
            cat("\n")

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            cat("  Iterations:", object@iterations, "\n\n")

            invisible(object)
          })

#' Calculate AIC for TS3fit object
#'
#' @param object TS3fit object
#' @param ... Additional arguments (not used)
#' @param k Penalty per parameter (default 2)
#'
#' @return AIC value
#' @export
setMethod("AIC", "TS3fit",
          function(object, ..., k = 2) {
            res <- object@residuals
            res <- res[is.finite(res)]
            n <- length(res)
            p <- length(object@coefficients)
            ll <- -n/2 * log(sum(res^2)/n) - n/2 * (1 + log(2*pi))
            -2 * ll + k * p
          })

#' Plot diagnostic plots for TS3fit object
#'
#' @param x TS3fit object
#' @param y Not used
#' @param which Set of plots to display (values 1-4)
#' @param ... Additional arguments
#'
#' @return Invisibly returns the input object
#' @export
setMethod("plot", signature(x = "TS3fit", y = "missing"),
          function(x, y, which = 1:4, ...) {
            res <- as.numeric(x@residuals)
            res <- res[is.finite(res)]

            which <- intersect(unique(which), 1:4)
            if (length(which) == 0) which <- 1:4

            old_par <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(old_par))
            graphics::par(mfrow = c(2, 2))

            for (idx in which) {
              switch(idx,
                     {
                       graphics::plot(seq_along(res), res, type = "l",
                                      main = "Residuals",
                                      xlab = "Index", ylab = "Residual", ...)
                       graphics::abline(h = 0, col = "red", lty = 2)
                     },
                     {
                       stats::qqnorm(res, main = "Normal Q-Q", ...)
                       stats::qqline(res, col = "red", lty = 2)
                     },
                     {
                       acf_vals <- stats::acf(res, plot = FALSE)
                       graphics::plot(acf_vals, main = "ACF of Residuals", ...)
                     },
                     {
                       graphics::hist(res, main = "Residual Histogram",
                                      xlab = "Residuals", breaks = "FD", ...)
                     })
            }

            invisible(x)
          })

#' Predict method for TS3fit objects
#'
#' @param object TS3fit object
#' @param n.ahead Integer: number of steps ahead to forecast
#' @param ... Additional arguments (not used)
#'
#' @return Numeric vector of predicted values
#' @export
setMethod("predict", "TS3fit",
          function(object, n.ahead = 1, ...) {
            if (object@model_type == "ar") {
              ar_order <- object@order$ar
              ar_coef  <- object@coefficients[1:ar_order]
              intercept <- object@intercept
              x <- object@original_series
              n <- length(x)
              x_centered <- if (intercept != 0) x - intercept else x

              preds <- numeric(n.ahead)
              history <- x_centered
              for (h in seq_len(n.ahead)) {
                val <- sum(ar_coef * rev(tail(history, ar_order)))
                preds[h] <- val + intercept
                history <- c(history, val)
              }
              return(preds)
            }

            # For MA/ARMA/ARIMA: use stats::arima with fixed params
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d

            arima_order <- c(ar_order, d, ma_order)
            fixed_params <- object@coefficients

            tryCatch({
              refit <- stats::arima(object@original_series,
                                    order = arima_order,
                                    fixed = fixed_params,
                                    include.mean = (object@intercept != 0))
              stats::predict(refit, n.ahead = n.ahead)$pred
            }, error = function(e) {
              warning("Prediction fallback: returning mean forecast. ", e$message)
              rep(mean(object@original_series), n.ahead)
            })
          })
