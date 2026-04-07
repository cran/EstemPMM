#' Calculate SARIMAX Residuals
#'
#' @param theta Combined vector of parameters (AR, MA, SAR, SMA, Intercept, Regressors)
#' @param y Time series data
#' @param xreg Exogenous regressors (optional)
#' @param order ARIMA order c(p, d, q)
#' @param seasonal Seasonal order c(P, D, Q)
#' @param period Seasonal period
#' @param include.mean Boolean, whether to include mean/intercept
#'
#' @return Vector of residuals
get_sarimax_residuals <- function(theta, y, xreg = NULL, order = c(0, 0, 0),
                                  seasonal = list(order = c(0, 0, 0), period = NA),
                                  include.mean = TRUE) {
    p <- order[1]
    d <- order[2]
    q <- order[3]
    P <- seasonal$order[1]
    D <- seasonal$order[2]
    Q <- seasonal$order[3]
    S <- seasonal$period

    # Parse theta
    idx <- 1
    ar_coefs <- if (p > 0) theta[idx:(idx + p - 1)] else numeric(0)
    idx <- idx + p
    ma_coefs <- if (q > 0) theta[idx:(idx + q - 1)] else numeric(0)
    idx <- idx + q
    sar_coefs <- if (P > 0) theta[idx:(idx + P - 1)] else numeric(0)
    idx <- idx + P
    sma_coefs <- if (Q > 0) theta[idx:(idx + Q - 1)] else numeric(0)
    idx <- idx + Q

    intercept <- if (include.mean) theta[idx] else 0
    if (include.mean) idx <- idx + 1

    xreg_coefs <- if (!is.null(xreg)) theta[idx:length(theta)] else numeric(0)

    # Construct model list for arima()
    # We use stats::arima to filter and get residuals.
    # Note: arima() definition of MA signs might differ from some texts, but it's consistent within R.

    fixed_params <- numeric(0)
    if (p > 0) fixed_params <- c(fixed_params, ar_coefs)
    if (q > 0) fixed_params <- c(fixed_params, ma_coefs)
    if (P > 0) fixed_params <- c(fixed_params, sar_coefs)
    if (Q > 0) fixed_params <- c(fixed_params, sma_coefs)
    if (include.mean) fixed_params <- c(fixed_params, intercept)
    if (!is.null(xreg)) fixed_params <- c(fixed_params, xreg_coefs)

    # arima() fails if parameters are non-stationary/non-invertible.
    # We return a penalty or large residuals if so.
    # However, for residuals calculation, we might want to use 'filter' directly if arima throws error.
    # Let's try using arima with fixed parameters first.

    tryCatch(
        {
            fit <- arima(y,
                order = order, seasonal = seasonal, xreg = xreg,
                include.mean = include.mean, fixed = fixed_params, transform.pars = FALSE
            )
            return(as.vector(residuals(fit)))
        },
        error = function(e) {
            # Fallback: Return large residuals to discourage optimizer
            return(rep(1e5, length(y)))
        }
    )
}

#' Calculate SARIMAX Jacobian (Numerical)
#'
#' @param theta Parameters
#' @param ... Arguments passed to get_sarimax_residuals
#'
#' @return Jacobian matrix (n x p)
get_sarimax_jacobian <- function(theta, ...) {
    # Use numDeriv::jacobian
    # The function passed to jacobian should return the vector of residuals
    func <- function(th) {
        get_sarimax_residuals(th, ...)
    }

    # Calculate Jacobian: d(residuals)/d(theta)
    # Note: In our unified_pmm2, we expect J = d(f)/d(theta) = -d(residuals)/d(theta)
    # So we negate the result from numDeriv
    J_res <- numDeriv::jacobian(func, theta)
    return(-J_res)
}
