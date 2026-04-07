# pmm3_ts_classes.R - S4 class hierarchy for PMM3 time series models

#' Base S4 class for PMM3 time series model results
#'
#' Stores results from time series estimation using PMM3 (S=3).
#' Designed for symmetric platykurtic innovations. Does NOT inherit
#' from BasePMM2 or TS2fit.
#'
#' @slot coefficients numeric vector of estimated parameters
#' @slot residuals numeric vector of final residuals/innovations
#' @slot m2 numeric second central moment of initial residuals
#' @slot m4 numeric fourth central moment of initial residuals
#' @slot m6 numeric sixth central moment of initial residuals
#' @slot gamma4 numeric excess kurtosis coefficient
#' @slot gamma6 numeric sixth-order cumulant coefficient
#' @slot g_coefficient numeric theoretical variance reduction factor g3
#' @slot kappa numeric moment ratio used in NR solver
#' @slot convergence logical whether algorithm converged
#' @slot iterations numeric number of iterations performed
#' @slot call original function call
#' @slot model_type character string indicating model type
#' @slot intercept numeric intercept value
#' @slot original_series numeric vector of original time series
#' @slot order list of order parameters
#'
#' @exportClass TS3fit
setClass("TS3fit",
         slots = c(coefficients    = "numeric",
                   residuals       = "numeric",
                   m2              = "numeric",
                   m4              = "numeric",
                   m6              = "numeric",
                   gamma4          = "numeric",
                   gamma6          = "numeric",
                   g_coefficient   = "numeric",
                   kappa           = "numeric",
                   convergence     = "logical",
                   iterations      = "numeric",
                   call            = "call",
                   model_type      = "character",
                   intercept       = "numeric",
                   original_series = "numeric",
                   order           = "list"))

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
