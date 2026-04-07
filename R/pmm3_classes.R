# pmm3_classes.R - S4 class definition for PMM3 models

#' S4 class for storing PMM3 regression model results
#'
#' PMM3 (S=3) is designed for symmetric platykurtic error distributions.
#' This class is fully standalone and does NOT inherit from BasePMM2.
#'
#' @slot coefficients numeric vector of estimated parameters
#' @slot residuals numeric vector of final residuals
#' @slot m2 numeric second central moment of initial residuals
#' @slot m4 numeric fourth central moment of initial residuals
#' @slot m6 numeric sixth central moment of initial residuals
#' @slot gamma4 numeric excess kurtosis coefficient
#' @slot gamma6 numeric sixth-order cumulant coefficient
#' @slot g_coefficient numeric theoretical variance reduction factor g3
#' @slot kappa numeric moment ratio used in NR solver
#' @slot convergence logical indicating whether algorithm converged
#' @slot iterations numeric number of iterations performed
#' @slot call original function call
#'
#' @exportClass PMM3fit
setClass("PMM3fit",
         slots = c(coefficients  = "numeric",
                   residuals     = "numeric",
                   m2            = "numeric",
                   m4            = "numeric",
                   m6            = "numeric",
                   gamma4        = "numeric",
                   gamma6        = "numeric",
                   g_coefficient = "numeric",
                   kappa         = "numeric",
                   convergence   = "logical",
                   iterations    = "numeric",
                   call          = "call"))
