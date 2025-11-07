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
