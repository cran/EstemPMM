# pmm_base_classes.R - Virtual base classes unifying the PMM hierarchy.
#
# Class diagram (introduced in 0.4.0; see EstemPMM_JSS Figure 1):
#
#   PMMfit (virtual)
#   |   coefficients, residuals, convergence, iterations, call
#   |
#   +-- BasePMM2 (virtual)        +-- BasePMM3 (virtual)
#   |     + m2, m3, m4            |     + m2, m4, m6,
#   |                             |       gamma4, gamma6,
#   |                             |       g_coefficient, kappa
#   |                             |
#   +-- PMMtsfit (virtual)
#         + model_type, intercept, original_series, order
#
# Concrete classes (PMM2fit, PMM3fit, TS2fit, TS3fit and their TS subclasses)
# inherit from these virtuals via single or multiple S4 inheritance. The
# diamond TS2fit -> {BasePMM2, PMMtsfit} -> PMMfit resolves cleanly because
# the only shared ancestor (PMMfit) is itself virtual with non-conflicting slots.

#' Virtual root S4 class for PMM fit objects
#'
#' Common parent of every concrete fit class returned by the package
#' (PMM2fit, PMM3fit and the time-series subclasses). Holds the slots
#' that are meaningful regardless of polynomial order or model family.
#'
#' This class is virtual and cannot be instantiated directly; it exists
#' so that S4 methods (e.g. \code{show}) can be defined once and inherited
#' by every concrete subclass.
#'
#' @slot coefficients numeric vector of estimated parameters
#' @slot residuals    numeric vector of final residuals/innovations
#' @slot convergence  logical or integer indicating algorithm convergence
#' @slot iterations   numeric number of iterations performed
#' @slot call         original function call
#'
#' @exportClass PMMfit
setClass("PMMfit",
         contains = "VIRTUAL",
         slots = c(coefficients = "numeric",
                   residuals    = "numeric",
                   convergence  = "logical",
                   iterations   = "numeric",
                   call         = "call"))

#' Virtual S4 class for PMM time-series fit objects
#'
#' Common parent of every concrete time-series fit class (TS2fit, TS3fit
#' and their AR/MA/ARMA/ARIMA/seasonal subclasses). Provides the
#' time-series-specific slots so that methods such as \code{predict},
#' \code{plot}, and \code{forecast} can be defined once and inherited.
#'
#' This class is virtual and inherits from \code{\linkS4class{PMMfit}}.
#'
#' @slot model_type      character string identifying the model family
#'   (one of \code{"ar"}, \code{"ma"}, \code{"arma"}, \code{"arima"},
#'   \code{"sar"}, \code{"sma"}, \code{"sarma"}, \code{"sarima"})
#' @slot intercept       numeric scalar intercept/mean term
#' @slot original_series numeric vector of the original time series
#' @slot order           list of order parameters (entries depend on
#'   \code{model_type}; see the relevant fitting function)
#'
#' @exportClass PMMtsfit
setClass("PMMtsfit",
         contains = c("PMMfit", "VIRTUAL"),
         slots = c(model_type      = "character",
                   intercept       = "numeric",
                   original_series = "numeric",
                   order           = "list"))
