# pmm_show_methods.R - Concise show/print methods for PMM fit objects.
#
# Replaces the default S4 slot-dump display with an lm()/arima()-style
# header so that typing `> fit` returns a human-readable summary.
# Detailed diagnostics remain available through summary().

#' @importFrom methods show
NULL

#' Show method for PMM2fit objects
#'
#' Prints a concise lm-style display: the original call, coefficients,
#' and a one-line algorithm footer. Use \code{\link{summary}} for full
#' diagnostics including moments, the g2 efficiency factor, and
#' bootstrap inference.
#'
#' @param object object of class \code{PMM2fit}
#'
#' @return The object (invisibly).
#'
#' @export
setMethod("show", "PMM2fit",
          function(object) {
            digits <- max(3L, getOption("digits") - 3L)
            cat("\nCall:\n",
                paste(deparse(object@call), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")
            if (length(object@coefficients)) {
              cat("Coefficients:\n")
              print.default(format(object@coefficients, digits = digits),
                            print.gap = 2L, quote = FALSE)
            } else {
              cat("No coefficients\n")
            }
            conv_str <- if (isTRUE(object@convergence)) "converged" else "did NOT converge"
            cat("\nPMM2 algorithm: ", conv_str,
                " in ", object@iterations, " iteration",
                if (object@iterations == 1L) "" else "s", ".\n", sep = "")
            invisible(object)
          })

#' Show method for PMM3fit objects
#'
#' Prints a concise lm-style display: the original call, coefficients,
#' and a one-line algorithm footer. Use \code{\link{summary}} for full
#' diagnostics including moments, gamma4/gamma6 cumulants, and the g3
#' efficiency factor.
#'
#' @param object object of class \code{PMM3fit}
#'
#' @return The object (invisibly).
#'
#' @export
setMethod("show", "PMM3fit",
          function(object) {
            digits <- max(3L, getOption("digits") - 3L)
            cat("\nCall:\n",
                paste(deparse(object@call), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")
            if (length(object@coefficients)) {
              cat("Coefficients:\n")
              print.default(format(object@coefficients, digits = digits),
                            print.gap = 2L, quote = FALSE)
            } else {
              cat("No coefficients\n")
            }
            conv_str <- if (isTRUE(object@convergence)) "converged" else "did NOT converge"
            cat("\nPMM3 algorithm: ", conv_str,
                " in ", object@iterations, " iteration",
                if (object@iterations == 1L) "" else "s", ".\n", sep = "")
            invisible(object)
          })

#' Show method for TS2fit objects (and subclasses)
#'
#' Inherited by ARPMM2, MAPMM2, ARMAPMM2, ARIMAPMM2 and the seasonal
#' SARPMM2/SMAPMM2/SARMAPMM2/SARIMAPMM2 subclasses. Use
#' \code{\link{summary}} for full diagnostics.
#'
#' @param object object inheriting from \code{TS2fit}
#'
#' @return The object (invisibly).
#'
#' @export
setMethod("show", "TS2fit",
          function(object) {
            digits <- max(3L, getOption("digits") - 3L)
            cat("\nCall:\n",
                paste(deparse(object@call), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")
            model_label <- .pmm_ts_model_label(object@model_type, object@order)
            cat("Model: ", model_label, "\n", sep = "")
            if (length(object@coefficients)) {
              cat("\nCoefficients:\n")
              print.default(format(object@coefficients, digits = digits),
                            print.gap = 2L, quote = FALSE)
            }
            if (length(object@intercept) && object@intercept != 0) {
              cat("\nIntercept: ", format(object@intercept, digits = digits), "\n", sep = "")
            }
            conv_str <- if (isTRUE(object@convergence)) "converged" else "did NOT converge"
            cat("\nPMM2 algorithm: ", conv_str,
                " in ", object@iterations, " iteration",
                if (object@iterations == 1L) "" else "s", ".\n", sep = "")
            invisible(object)
          })

#' Show method for TS3fit objects (and subclasses)
#'
#' Inherited by ARPMM3, MAPMM3, ARMAPMM3, ARIMAPMM3. Use
#' \code{\link{summary}} for full diagnostics.
#'
#' @param object object inheriting from \code{TS3fit}
#'
#' @return The object (invisibly).
#'
#' @export
setMethod("show", "TS3fit",
          function(object) {
            digits <- max(3L, getOption("digits") - 3L)
            cat("\nCall:\n",
                paste(deparse(object@call), sep = "\n", collapse = "\n"),
                "\n\n", sep = "")
            model_label <- .pmm_ts_model_label(object@model_type, object@order)
            cat("Model: ", model_label, "\n", sep = "")
            if (length(object@coefficients)) {
              cat("\nCoefficients:\n")
              print.default(format(object@coefficients, digits = digits),
                            print.gap = 2L, quote = FALSE)
            }
            if (length(object@intercept) && object@intercept != 0) {
              cat("\nIntercept: ", format(object@intercept, digits = digits), "\n", sep = "")
            }
            conv_str <- if (isTRUE(object@convergence)) "converged" else "did NOT converge"
            cat("\nPMM3 algorithm: ", conv_str,
                " in ", object@iterations, " iteration",
                if (object@iterations == 1L) "" else "s", ".\n", sep = "")
            invisible(object)
          })

# Helper: format a TS model label like "ARIMA(2,1,1)" or
# "SARIMA(1,0,0)x(1,1,0)[12]" from the @model_type tag and @order list.
.pmm_ts_model_label <- function(model_type, order) {
  ord <- function(name, default = 0L) {
    val <- order[[name]]
    if (is.null(val)) default else as.integer(val)
  }
  type <- tolower(model_type)
  switch(type,
    "ar"     = sprintf("AR(%d)",     ord("ar")),
    "ma"     = sprintf("MA(%d)",     ord("ma")),
    "arma"   = sprintf("ARMA(%d,%d)", ord("ar"), ord("ma")),
    "arima"  = sprintf("ARIMA(%d,%d,%d)", ord("ar"), ord("d"), ord("ma")),
    "sar"    = sprintf("SAR(%d)x(%d)[%d]",
                       ord("ar"), ord("sar"), ord("period")),
    "sma"    = sprintf("SMA(%d)x(%d)[%d]",
                       ord("ma"), ord("sma"), ord("period")),
    "sarma"  = sprintf("SARMA(%d,%d)x(%d,%d)[%d]",
                       ord("ar"), ord("ma"), ord("sar"), ord("sma"), ord("period")),
    "sarima" = sprintf("SARIMA(%d,%d,%d)x(%d,%d,%d)[%d]",
                       ord("ar"), ord("d"), ord("ma"),
                       ord("sar"), ord("D"), ord("sma"), ord("period")),
    toupper(type)
  )
}
