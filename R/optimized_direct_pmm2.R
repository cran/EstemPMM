source("R/sarimax_wrapper.R")

#' @title Optimized Direct Nonlinear PMM2 (DEPRECATED)
#' @description
#' \strong{DEPRECATED:} This function has been deprecated and should not be used.
#'
#' \strong{Reason for Deprecation:}
#' Monte Carlo experiments (R=50, n=200) showed this direct optimization approach
#' performs catastrophically worse than classical methods:
#' \itemize{
#'   \item AR(1): 17.0x worse MSE than MLE (3.67 vs 0.22)
#'   \item MA(1): Failed to converge in most replications
#'   \item SARIMA: Severe numerical instability
#' }
#'
#' \strong{Recommended Alternatives:}
#' Use the Unified PMM2 framework instead, which showed 3-23\% MSE improvements:
#' \itemize{
#'   \item \code{\link{ar_pmm2}} with \code{pmm2_variant = "unified_global"} (default)
#'   \item \code{\link{arima_pmm2}} with \code{pmm2_variant = "unified_iterative"}
#'   \item \code{\link{ma_pmm2}} with \code{pmm2_variant = "linearized"}
#' }
#'
#' @details
#' This function attempted to directly maximize the PMM2 objective function:
#' \deqn{J(\theta) = \gamma_3^2 / (2 + \gamma_4)}
#'
#' However, this objective surface is highly non-convex with many local optima
#' and flat regions, making direct optimization unreliable. The Unified PMM2
#' approach (global/iterative corrections from classical estimates) proved far
#' more stable and accurate.
#'
#' \strong{Technical Issues:}
#' \itemize{
#'   \item Objective function has flat regions near optimal values
#'   \item Gradient becomes numerically unstable when \eqn{\gamma_4 \approx -2}
#'   \item Sensitive to initialization (even from MLE estimates)
#'   \item No guarantee of convergence to global optimum
#' }
#'
#' @param theta_init Initial parameters (not used in deprecated version)
#' @param y Time series (not used in deprecated version)
#' @param order ARIMA order (not used in deprecated version)
#' @param seasonal Seasonal order (not used in deprecated version)
#'
#' @return This function will throw an error directing users to unified methods
#'
#' @references
#' Monte Carlo comparison results documented in NEWS.md (Version 0.2.0)
#'
#' @keywords internal
#'
#' @seealso \code{\link{ar_pmm2}}, \code{\link{ma_pmm2}}, \code{\link{arima_pmm2}}
#'
#' Optimized Direct Nonlinear PMM2
#'
#' Uses analytical gradients to maximize the PMM2 objective function.
#' J(theta) = gamma3^2 / (2 + gamma4)
#'
#' @param theta_init Initial parameters
#' @param y Time series
#' @param order ARIMA order
#' @param seasonal Seasonal order
#' @return Optimization result
optimized_direct_pmm2 <- function(theta_init, y, order = c(0, 0, 0), seasonal = list(order = c(0, 0, 0), period = NA)) {
    # DEPRECATED: This function fails catastrophically in practice (17x worse MSE than MLE)
    # Use unified PMM2 variants instead: ar_pmm2(), ma_pmm2(), arima_pmm2() with pmm2_variant parameter
    .Deprecated(
        msg = paste(
            "optimized_direct_pmm2() is deprecated due to poor performance (17x worse MSE than MLE).",
            "Use unified PMM2 variants instead:",
            "  - ar_pmm2(x, order=..., pmm2_variant='unified_global')",
            "  - ma_pmm2(x, order=..., pmm2_variant='linearized')",
            "  - arima_pmm2(x, order=..., pmm2_variant='unified_iterative')",
            "See NEWS.md (Version 0.2.0) for Monte Carlo comparison results.",
            sep = "\n"
        )
    )

    stop(
        "optimized_direct_pmm2() has been deprecated. ",
        "Use ar_pmm2(), ma_pmm2(), or arima_pmm2() with pmm2_variant parameter instead. ",
        "See ?ar_pmm2 for details."
    )

    # Original implementation commented out (preserved for historical reference)
    # Wrapper for objective and gradient
    # fn_obj_grad <- function(theta) {
    #     # 1. Get Residuals and Jacobian
    #     # Note: get_sarimax_residuals returns e_t
    #     # get_sarimax_jacobian returns J = -d(e)/d(theta) (as per unified definition)
    #     # But for gradient of moments we need d(e)/d(theta) directly.
    #     # Let's use the wrapper but be careful with signs.
    #
    #     e <- get_sarimax_residuals(theta, y, order = order, seasonal = seasonal)
    #
    #     # Check for bad parameters (arima returns large residuals)
    #     if (any(abs(e) > 1e4)) {
    #         return(list(value = -1e10, gradient = rep(0, length(theta))))
    #     }
    #
    #     # Jacobian from wrapper is J_unified = - d(e)/d(theta)
    #     # So d(e)/d(theta) = - J_unified
    #     J_unified <- get_sarimax_jacobian(theta, y, order = order, seasonal = seasonal)
    #     d_e <- -J_unified
    #
    #     n <- length(e)
    #
    #     # 2. Compute Moments
    #     m2 <- mean(e^2)
    #     m3 <- mean(e^3)
    #     m4 <- mean(e^4)
    #
    #     if (m2 < 1e-8) {
    #         return(list(value = -1e10, gradient = rep(0, length(theta))))
    #     }
    #
    #     gamma3 <- m3 / m2^(1.5)
    #     gamma4 <- m4 / m2^2 - 3
    #
    #     denom <- 2 + gamma4
    #     if (denom < 0.01) denom <- 0.01 # Regularization
    #
    #     obj <- gamma3^2 / denom
    #
    #     # 3. Compute Gradients of Moments
    #     # d(m_k) / d(theta) = 1/n * sum( k * e^(k-1) * d_e/d_theta )
    #     # d_e is (n x p) matrix
    #
    #     # Helper to compute grad of mean(e^k)
    #     grad_moment <- function(k) {
    #         # vec = k * e^(k-1)
    #         # result = 1/n * t(d_e) %*% vec
    #         vec <- k * e^(k - 1)
    #         as.vector(t(d_e) %*% vec) / n
    #     }
    #
    #     dm2 <- grad_moment(2)
    #     dm3 <- grad_moment(3)
    #     dm4 <- grad_moment(4)
    #
    #     # 4. Gradients of Gammas
    #     # gamma3 = m3 * m2^(-1.5)
    #     # d(gamma3) = dm3 * m2^(-1.5) + m3 * (-1.5) * m2^(-2.5) * dm2
    #     d_gamma3 <- dm3 * m2^(-1.5) - 1.5 * m3 * m2^(-2.5) * dm2
    #
    #     # gamma4 = m4 * m2^(-2) - 3
    #     # d(gamma4) = dm4 * m2^(-2) + m4 * (-2) * m2^(-3) * dm2
    #     d_gamma4 <- dm4 * m2^(-2) - 2 * m4 * m2^(-3) * dm2
    #
    #     # 5. Gradient of Objective
    #     # J = g3^2 / (2 + g4)
    #     # dJ = [ 2*g3*dg3 * (2+g4) - g3^2 * dg4 ] / (2+g4)^2
    #
    #     grad_obj <- (2 * gamma3 * d_gamma3 * denom - gamma3^2 * d_gamma4) / denom^2
    #
    #     # Return list for optim (if using check_grad) or just value/grad depending on method
    #     # optim expects fn to return value, gr to return gradient
    #     # We will separate them below
    #     attr(obj, "gradient") <- grad_obj
    #     obj
    # }
    #
    # fn_val <- function(theta) {
    #     res <- fn_obj_grad(theta)
    #     as.numeric(res)
    # }
    #
    # fn_grad <- function(theta) {
    #     res <- fn_obj_grad(theta)
    #     attr(res, "gradient")
    # }
    #
    # # Optimization
    # # Use L-BFGS-B with bounds to prevent explosion
    # lower_bounds <- rep(-0.99, length(theta_init))
    # upper_bounds <- rep(0.99, length(theta_init))
    #
    # # Adjust bounds for intercept (no bounds)
    # if ("intercept" %in% names(theta_init)) {
    #     lower_bounds[names(theta_init) == "intercept"] <- -Inf
    #     upper_bounds[names(theta_init) == "intercept"] <- Inf
    # }
    #
    # tryCatch(
    #     {
    #         res <- optim(theta_init, fn_val,
    #             gr = fn_grad, method = "L-BFGS-B",
    #             lower = lower_bounds, upper = upper_bounds,
    #             control = list(fnscale = -1, maxit = 100)
    #         ) # fnscale -1 for maximization
    #         list(coefficients = res$par, convergence = res$convergence, value = res$value)
    #     },
    #     error = function(e) {
    #         list(coefficients = theta_init, convergence = 1, value = NA, error = e$message)
    #     }
    # )
}

