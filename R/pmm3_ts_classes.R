# pmm3_ts_classes.R - S4 class hierarchy for PMM3 time-series models.

#' @include pmm_base_classes.R pmm3_classes.R
NULL

#' S4 class for PMM3 time-series fit results
#'
#' Common parent of every PMM3 time-series fit class (ARPMM3, MAPMM3,
#' ARMAPMM3, ARIMAPMM3). Inherits the PMM3-specific moments and cumulant
#' coefficients from \code{\linkS4class{BasePMM3}} and the time-series
#' slots from \code{\linkS4class{PMMtsfit}}. The diamond inheritance
#' resolves at \code{\linkS4class{PMMfit}}, which both branches share
#' as their virtual root.
#'
#' @exportClass TS3fit
setClass("TS3fit",
         contains = c("BasePMM3", "PMMtsfit"))

#' S4 class for PMM3 AR model results
#' @exportClass ARPMM3
setClass("ARPMM3", contains = "TS3fit")

#' S4 class for PMM3 MA model results
#' @exportClass MAPMM3
setClass("MAPMM3", contains = "TS3fit")

#' S4 class for PMM3 ARMA model results
#' @exportClass ARMAPMM3
setClass("ARMAPMM3", contains = "TS3fit")

#' S4 class for PMM3 ARIMA model results
#' @exportClass ARIMAPMM3
setClass("ARIMAPMM3", contains = "TS3fit")
