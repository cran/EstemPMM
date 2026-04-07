# pmm2_classes.R - Class hierarchy for PMM2 models

#' Base S4 class for storing PMM2 model results
#'
#' @slot coefficients numeric vector of estimated parameters
#' @slot residuals numeric vector of final residuals
#' @slot m2 numeric second central moment of initial residuals
#' @slot m3 numeric third central moment of initial residuals
#' @slot m4 numeric fourth central moment of initial residuals
#' @slot convergence logical or integer code indicating whether algorithm converged
#' @slot iterations numeric number of iterations performed
#' @slot call original function call
#'
#' @exportClass BasePMM2
setClass("BasePMM2",
         slots = c(coefficients = "numeric",
                   residuals    = "numeric",
                   m2           = "numeric",
                   m3           = "numeric",
                   m4           = "numeric",
                   convergence  = "logical",
                   iterations   = "numeric",
                   call         = "call"))

#' S4 class for storing PMM2 regression model results
#'
#' @slot coefficients numeric vector of estimated parameters
#' @slot residuals numeric vector of final residuals
#' @slot m2 numeric second central moment of initial residuals
#' @slot m3 numeric third central moment of initial residuals
#' @slot m4 numeric fourth central moment of initial residuals
#' @slot convergence logical or integer code indicating whether algorithm converged
#' @slot iterations numeric number of iterations performed
#' @slot call original function call
#'
#' @exportClass PMM2fit
setClass("PMM2fit",
         contains = "BasePMM2")

#' Base S4 class for storing PMM2 time series model results
#'
#' @slot coefficients numeric vector of estimated parameters
#' @slot residuals numeric vector of final residuals
#' @slot m2 numeric second central moment of initial residuals
#' @slot m3 numeric third central moment of initial residuals
#' @slot m4 numeric fourth central moment of initial residuals
#' @slot convergence logical or integer code indicating whether algorithm converged
#' @slot iterations numeric number of iterations performed
#' @slot call original function call
#' @slot model_type character string indicating model type
#' @slot intercept numeric value of intercept
#' @slot original_series numeric vector of original time series
#' @slot order list of order parameters
#'
#' @exportClass TS2fit
setClass("TS2fit",
         contains = "BasePMM2",
         slots = c(model_type      = "character",
                   intercept       = "numeric",
                   original_series = "numeric",
                   order           = "list"))

#' S4 class for storing PMM2 AR model results
#'
#' @exportClass ARPMM2
setClass("ARPMM2", contains = "TS2fit")

#' S4 class for storing PMM2 MA model results
#'
#' @exportClass MAPMM2
setClass("MAPMM2", contains = "TS2fit")

#' S4 class for storing PMM2 ARMA model results
#'
#' @exportClass ARMAPMM2
setClass("ARMAPMM2", contains = "TS2fit")

#' S4 class for storing PMM2 ARIMA model results
#'
#' @exportClass ARIMAPMM2
setClass("ARIMAPMM2", contains = "TS2fit")

#' S4 class for Seasonal AR model results with PMM2
#'
#' This class stores the results of fitting a Seasonal Autoregressive (SAR)
#' model using the PMM2 method. It extends the TS2fit class with additional
#' slots specific to seasonal models.
#'
#' @slot coefficients Numeric vector of estimated parameters (AR and SAR coefficients)
#' @slot residuals Numeric vector of residuals/innovations
#' @slot m2 Second central moment of residuals
#' @slot m3 Third central moment of residuals
#' @slot m4 Fourth central moment of residuals
#' @slot convergence Logical, whether PMM2 algorithm converged
#' @slot iterations Integer, number of iterations performed
#' @slot call Original function call
#' @slot model_type Character, model type identifier ("sar")
#' @slot intercept Numeric, intercept/mean term
#' @slot original_series Numeric vector, original time series data
#' @slot order List with model specification: list(ar, sar, period)
#'   \itemize{
#'     \item ar: Non-seasonal AR order (p)
#'     \item sar: Seasonal AR order (P)
#'     \item period: Seasonal period (s)
#'   }
#'
#' @details
#' The SARPMM2 class represents fitted SAR models of the form
#' \deqn{y_t = \phi_1 y_{t-1} + \dots + \phi_p y_{t-p} +
#'             \Phi_1 y_{t-s} + \dots + \Phi_P y_{t-Ps} + \epsilon_t.}
#'
#' Where:
#'   - p is the non-seasonal AR order
#'   - P is the seasonal AR order
#'   - s is the seasonal period
#'
#' @seealso \code{\link{sar_pmm2}} for fitting SAR models
#'
#' @exportClass SARPMM2
setClass("SARPMM2",
         slots = c(
           coefficients = "numeric",
           residuals = "numeric",
           m2 = "numeric",
           m3 = "numeric",
           m4 = "numeric",
           convergence = "logical",
           iterations = "numeric",
           call = "call",
           model_type = "character",
           intercept = "numeric",
           original_series = "numeric",
           order = "list"
         ),
         contains = "TS2fit")

#' S4 class for Seasonal MA PMM2 results
#'
#' This class stores the results of fitting a Seasonal Moving Average (SMA)
#' model using the Polynomial Maximization Method (PMM2).
#'
#' @slot coefficients Estimated seasonal MA coefficients (Theta_1, Theta_2, ..., Theta_Q)
#' @slot innovations Estimated innovations (residuals epsilon_t)
#' @slot m2 Second central moment (variance) of innovations
#' @slot m3 Third central moment (skewness indicator) of innovations
#' @slot m4 Fourth central moment (kurtosis indicator) of innovations
#' @slot convergence Logical indicating whether PMM2 algorithm converged
#' @slot iterations Number of iterations required for convergence
#' @slot call The function call that created this object
#' @slot model_type Character string "sma"
#' @slot intercept Model intercept (mean)
#' @slot original_series Original time series data
#' @slot order List with Q (seasonal MA order) and s (seasonal period)
#'
#' @details
#' The SMA(Q)_s model is expressed as
#' \deqn{y_t = \mu + \epsilon_t + \Theta_1 \epsilon_{t-s} + \Theta_2 \epsilon_{t-2s} + \dots + \Theta_Q \epsilon_{t-Qs}.}
#'
#' Where:
#'   - Q is the seasonal MA order
#'   - s is the seasonal period
#'   - epsilon_t are innovations
#'
#' @seealso \code{\link{sma_pmm2}} for fitting SMA models
#'
#' @exportClass SMAPMM2
setClass("SMAPMM2",
         slots = c(
           coefficients = "numeric",
           innovations = "numeric",
           m2 = "numeric",
           m3 = "numeric",
           m4 = "numeric",
           convergence = "logical",
           iterations = "numeric",
           call = "call",
           model_type = "character",
           intercept = "numeric",
           original_series = "numeric",
           order = "list"
         ),
         contains = "TS2fit")

#' S4 class for Seasonal ARMA model results with PMM2
#'
#' This class stores the results of fitting a Seasonal Autoregressive Moving Average
#' (SARMA) model using the PMM2 method. It combines both seasonal AR and seasonal MA
#' components.
#'
#' @slot coefficients Numeric vector of estimated parameters (SAR and SMA coefficients)
#' @slot residuals Numeric vector of residuals/innovations
#' @slot m2 Second central moment of residuals
#' @slot m3 Third central moment of residuals
#' @slot m4 Fourth central moment of residuals
#' @slot convergence Logical, whether PMM2 algorithm converged
#' @slot iterations Integer, number of iterations performed
#' @slot call Original function call
#' @slot model_type Character, model type identifier ("sarma")
#' @slot intercept Numeric, intercept/mean term
#' @slot original_series Numeric vector, original time series data
#' @slot order List with model specification: list(ar, sar, ma, sma, period)
#'   \itemize{
#'     \item ar: Non-seasonal AR order (p)
#'     \item sar: Seasonal AR order (P)
#'     \item ma: Non-seasonal MA order (q)
#'     \item sma: Seasonal MA order (Q)
#'     \item period: Seasonal period (s)
#'   }
#'
#' @details
#' The SARMAPMM2 class represents fitted SARMA models combining:
#' \itemize{
#'   \item AR(p): \eqn{\phi_1 y_{t-1} + \dots + \phi_p y_{t-p}}
#'   \item Seasonal AR component: \eqn{\Phi_1 y_{t-s} + \dots + \Phi_P y_{t-Ps}}
#'   \item MA(q): \eqn{\theta_1 \epsilon_{t-1} + \dots + \theta_q \epsilon_{t-q}}
#'   \item Seasonal MA component: \eqn{\Theta_1 \epsilon_{t-s} + \dots + \Theta_Q \epsilon_{t-Qs}}
#' }
#'
#' @seealso \code{\link{sarma_pmm2}} for fitting SARMA models
#'
#' @exportClass SARMAPMM2
setClass("SARMAPMM2",
         slots = c(
           coefficients = "numeric",
           residuals = "numeric",
           m2 = "numeric",
           m3 = "numeric",
           m4 = "numeric",
           convergence = "logical",
           iterations = "numeric",
           call = "call",
           model_type = "character",
           intercept = "numeric",
           original_series = "numeric",
           order = "list"
         ),
         contains = "TS2fit")

#' S4 class for Seasonal ARIMA model results with PMM2
#'
#' This class stores the results of fitting a Seasonal Autoregressive Integrated
#' Moving Average (SARIMA) model using the PMM2 method. It extends SARMA with
#' differencing operators.
#'
#' @slot coefficients Numeric vector of estimated parameters
#' @slot residuals Numeric vector of residuals/innovations
#' @slot m2 Second central moment of residuals
#' @slot m3 Third central moment of residuals
#' @slot m4 Fourth central moment of residuals
#' @slot convergence Logical, whether PMM2 algorithm converged
#' @slot iterations Integer, number of iterations performed
#' @slot call Original function call
#' @slot model_type Character, model type identifier ("sarima")
#' @slot intercept Numeric, intercept/mean term
#' @slot original_series Numeric vector, original time series data
#' @slot order List with model specification: list(ar, sar, ma, sma, d, D, period)
#'   \itemize{
#'     \item ar: Non-seasonal AR order (p)
#'     \item sar: Seasonal AR order (P)
#'     \item ma: Non-seasonal MA order (q)
#'     \item sma: Seasonal MA order (Q)
#'     \item d: Non-seasonal differencing order
#'     \item D: Seasonal differencing order
#'     \item period: Seasonal period (s)
#'   }
#'
#' @details
#' The SARIMAPMM2 class represents fitted SARIMA(p,d,q) x (P,D,Q)_s models:
#' \deqn{(1 - \phi_1 B - \dots - \phi_p B^p)(1 - \Phi_1 B^s - \dots - \Phi_P B^{Ps})(1 - B)^d (1 - B^s)^D y_t = (1 + \theta_1 B + \dots + \theta_q B^q)(1 + \Theta_1 B^s + \dots + \Theta_Q B^{Qs}) \epsilon_t.}
#'
#' Where B is the backshift operator.
#'
#' @seealso \code{\link{sarima_pmm2}} for fitting SARIMA models
#'
#' @exportClass SARIMAPMM2
setClass("SARIMAPMM2",
         slots = c(
           coefficients = "numeric",
           residuals = "numeric",
           m2 = "numeric",
           m3 = "numeric",
           m4 = "numeric",
           convergence = "logical",
           iterations = "numeric",
           call = "call",
           model_type = "character",
           intercept = "numeric",
           original_series = "numeric",
           order = "list"
         ),
         contains = "TS2fit")

# =============================================================================
# Methods
# =============================================================================

#' Generic summary method for PMM2fit objects
#'
#' @param object object of class "PMM2fit"
#' @param formula (optional) formula used for the model
#' @param data (optional) data used
#' @param B number of bootstrap replications for statistical inference
#' @param ... additional arguments (not used)
#'
#' @return Prints summary to console; returns object (invisibly).
#'
#' @export
setMethod("summary", "PMM2fit",
          function(object, formula=NULL, data=NULL, B=100, ...) {
            cat("PMM2 estimation results\n")
            if(!is.null(object@call)) {
              cat("Call:\n")
              print(object@call)
              cat("\n")
            }

            cat("Coefficients:\n")
            print(object@coefficients)

            cat("\nCentral moments of initial residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Theoretical characteristics of PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (expected ratio Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            if("iterations" %in% slotNames(object)) {
              cat("  Iterations:", object@iterations, "\n\n")
            } else {
              cat("\n")
            }

            # If user wants to see p-values, call pmm2_inference:
            if(!is.null(formula) && !is.null(data)) {
              cat("Approximate statistical inference via bootstrap (B=", B, "):\n", sep="")
              inf_tab <- pmm2_inference(object, formula, data, B=B)
              print(inf_tab)
            } else {
              cat("To view p-values, provide formula= and data=\n")
            }
            invisible(object)
          }
)

#' Generic summary method for TS2fit objects
#'
#' @param object object of class "TS2fit" or subclass
#' @param ... additional arguments (not used)
#'
#' @return Prints summary to console; returns object (invisibly).
#'
#' @export
setMethod("summary", "TS2fit",
          function(object, ...) {
            model_type <- object@model_type
            ar_order <- object@order$ar
            ma_order <- object@order$ma
            d <- object@order$d

            cat("PMM2 time series estimation results\n")

            # Print model type
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

            if(!is.null(object@call)) {
              cat("Call:\n")
              print(object@call)
              cat("\n")
            }

            # Print coefficients with separation into AR and MA parts
            cat("Coefficients:\n")
            if (ar_order > 0) {
              ar_coefs <- object@coefficients[1:ar_order]
              names(ar_coefs) <- paste0("ar", 1:ar_order)
              cat("AR: ")
              print(ar_coefs)
            }

            if (ma_order > 0) {
              ma_coefs <- object@coefficients[(ar_order+1):(ar_order+ma_order)]
              names(ma_coefs) <- paste0("ma", 1:ma_order)
              cat("MA: ")
              print(ma_coefs)
            }

            # Print intercept if present
            if (object@intercept != 0) {
              cat("Intercept: ", object@intercept, "\n")
            }

            cat("\nCentral moments of initial residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Theoretical characteristics of PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (expected ratio Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            if("iterations" %in% slotNames(object)) {
              cat("  Iterations:", object@iterations, "\n\n")
            } else {
              cat("\n")
            }

            invisible(object)
          }
)

#' Generic summary method for SARMAPMM2 objects
#'
#' @param object object of class "SARMAPMM2"
#' @param ... additional arguments (not used)
#'
#' @return Prints summary to console; returns object (invisibly).
#'
#' @export
setMethod("summary", "SARMAPMM2",
          function(object, ...) {
            p <- object@order$ar
            P <- object@order$sar
            q <- object@order$ma
            Q <- object@order$sma
            s <- object@order$period

            cat("PMM2 Seasonal ARMA estimation results\n")
            cat("Model type: SARMA(", p, ",", q, ")x(", P, ",", Q, ")_", s, "\n", sep = "")

            if(!is.null(object@call)) {
              cat("Call:\n")
              print(object@call)
              cat("\n")
            }

            cat("Coefficients:\n")
            coef_idx <- 1

            if (p > 0) {
              ar_coefs <- object@coefficients[coef_idx:(coef_idx + p - 1)]
              names(ar_coefs) <- paste0("ar", 1:p)
              cat("AR: ")
              print(ar_coefs)
              coef_idx <- coef_idx + p
            }

            if (P > 0) {
              sar_coefs <- object@coefficients[coef_idx:(coef_idx + P - 1)]
              names(sar_coefs) <- paste0("sar", 1:P)
              cat("SAR: ")
              print(sar_coefs)
              coef_idx <- coef_idx + P
            }

            if (q > 0) {
              ma_coefs <- object@coefficients[coef_idx:(coef_idx + q - 1)]
              names(ma_coefs) <- paste0("ma", 1:q)
              cat("MA: ")
              print(ma_coefs)
              coef_idx <- coef_idx + q
            }

            if (Q > 0) {
              sma_coefs <- object@coefficients[coef_idx:(coef_idx + Q - 1)]
              names(sma_coefs) <- paste0("sma", 1:Q)
              cat("SMA: ")
              print(sma_coefs)
            }

            if (object@intercept != 0) {
              cat("Intercept: ", object@intercept, "\n")
            }

            cat("\nCentral moments of residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Theoretical characteristics of PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (expected ratio Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            cat("  Iterations:", object@iterations, "\n\n")

            invisible(object)
          }
)

#' Generic summary method for SARIMAPMM2 objects
#'
#' @param object object of class "SARIMAPMM2"
#' @param ... additional arguments (not used)
#'
#' @return Prints summary to console; returns object (invisibly).
#'
#' @export
setMethod("summary", "SARIMAPMM2",
          function(object, ...) {
            p <- object@order$ar
            P <- object@order$sar
            q <- object@order$ma
            Q <- object@order$sma
            d <- object@order$d
            D <- object@order$D
            s <- object@order$period

            cat("PMM2 Seasonal ARIMA estimation results\n")
            cat("Model type: SARIMA(", p, ",", d, ",", q, ")x(", P, ",", D, ",", Q, ")_", s, "\n", sep = "")

            if(!is.null(object@call)) {
              cat("Call:\n")
              print(object@call)
              cat("\n")
            }

            cat("Coefficients:\n")
            coef_idx <- 1

            if (p > 0) {
              ar_coefs <- object@coefficients[coef_idx:(coef_idx + p - 1)]
              names(ar_coefs) <- paste0("ar", 1:p)
              cat("AR: ")
              print(ar_coefs)
              coef_idx <- coef_idx + p
            }

            if (P > 0) {
              sar_coefs <- object@coefficients[coef_idx:(coef_idx + P - 1)]
              names(sar_coefs) <- paste0("sar", 1:P)
              cat("SAR: ")
              print(sar_coefs)
              coef_idx <- coef_idx + P
            }

            if (q > 0) {
              ma_coefs <- object@coefficients[coef_idx:(coef_idx + q - 1)]
              names(ma_coefs) <- paste0("ma", 1:q)
              cat("MA: ")
              print(ma_coefs)
              coef_idx <- coef_idx + q
            }

            if (Q > 0) {
              sma_coefs <- object@coefficients[coef_idx:(coef_idx + Q - 1)]
              names(sma_coefs) <- paste0("sma", 1:Q)
              cat("SMA: ")
              print(sma_coefs)
            }

            if (object@intercept != 0) {
              cat("Intercept: ", object@intercept, "\n")
            }

            cat("\nCentral moments of residuals:\n")
            cat("  m2 =", object@m2, "\n")
            cat("  m3 =", object@m3, "\n")
            cat("  m4 =", object@m4, "\n\n")

            vf <- pmm2_variance_factor(object@m2, object@m3, object@m4)
            if(!is.na(vf$g)) {
              cat("Theoretical characteristics of PMM2 (S = 2):\n")
              cat("  c3 =", vf$c3, "\n")
              cat("  c4 =", vf$c4, "\n")
              cat("  g  =", vf$g, " (expected ratio Var[PMM2]/Var[OLS])\n\n")
            }

            cat("Algorithm information:\n")
            cat("  Convergence status:", object@convergence, "\n")
            cat("  Iterations:", object@iterations, "\n\n")

            invisible(object)
          }
)
