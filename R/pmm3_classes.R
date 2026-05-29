# pmm3_classes.R - S4 class definitions for PMM3 models.

#' @include pmm_base_classes.R
NULL

#' Virtual S4 class for the PMM3 model family
#'
#' Virtual intermediate that inherits common slots from
#' \code{\linkS4class{PMMfit}} and adds the PMM3-specific moments and
#' cumulant coefficients used for symmetric platykurtic errors.
#' Concrete subclasses are \code{\linkS4class{PMM3fit}} (regression)
#' and the time-series classes under \code{\linkS4class{TS3fit}}.
#'
#' @slot m2            numeric second central moment of initial residuals
#' @slot m4            numeric fourth central moment of initial residuals
#' @slot m6            numeric sixth central moment of initial residuals
#' @slot gamma4        numeric excess kurtosis coefficient
#' @slot gamma6        numeric sixth-order cumulant coefficient
#' @slot g_coefficient numeric theoretical variance reduction factor g3
#' @slot kappa         numeric moment ratio used by the Newton-Raphson solver
#'
#' @exportClass BasePMM3
setClass("BasePMM3",
         contains = c("PMMfit", "VIRTUAL"),
         slots = c(m2            = "numeric",
                   m4            = "numeric",
                   m6            = "numeric",
                   gamma4        = "numeric",
                   gamma6        = "numeric",
                   g_coefficient = "numeric",
                   kappa         = "numeric"))

#' S4 class for PMM3 regression fit results
#'
#' Concrete class returned by \code{\link{lm_pmm3}}. Inherits all slots
#' from \code{\linkS4class{BasePMM3}}, which in turn inherits common
#' slots from \code{\linkS4class{PMMfit}}.
#'
#' @exportClass PMM3fit
setClass("PMM3fit",
         contains = "BasePMM3")
