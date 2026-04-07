# pmm3_dispatch.R - Automatic PMM method selection

#' Compute residual cumulants up to 6th order
#' @param residuals Numeric vector of residuals
#' @return Named list: m2, m3, m4, m6, gamma3, gamma4, gamma6
#' @keywords internal
.residual_cumulants <- function(residuals) {
  r  <- residuals - mean(residuals)
  m2 <- mean(r^2)
  m3 <- mean(r^3)
  m4 <- mean(r^4)
  m6 <- mean(r^6)

  gamma3 <- m3 / m2^(3/2)
  gamma4 <- m4 / m2^2 - 3
  gamma6 <- m6 / m2^3 - 15 * (m4 / m2^2) + 30

  list(m2 = m2, m3 = m3, m4 = m4, m6 = m6,
       gamma3 = gamma3, gamma4 = gamma4, gamma6 = gamma6)
}

#' PMM2 theoretical efficiency factor
#' @param gamma3 Skewness coefficient
#' @param gamma4 Excess kurtosis coefficient
#' @return Numeric scalar: g2 factor
#' @keywords internal
.g2_factor <- function(gamma3, gamma4) {
  denom <- 2 + gamma4
  if (is.na(denom) || denom <= 0) return(1.0)
  1 - gamma3^2 / denom
}

#' PMM3 theoretical efficiency factor
#' @param gamma4 Excess kurtosis coefficient
#' @param gamma6 Sixth-order cumulant coefficient
#' @return Numeric scalar: g3 factor
#' @keywords internal
.g3_factor <- function(gamma4, gamma6) {
  denom <- 6 + 9 * gamma4 + gamma6
  if (is.na(denom) || denom <= 0) return(1.0)
  1 - gamma4^2 / denom
}

#' Automatic PMM method selection
#'
#' Analyses OLS residual cumulants to recommend the best estimation method:
#' OLS (Gaussian errors), PMM2 (asymmetric errors), or PMM3 (symmetric
#' platykurtic errors).
#'
#' @param residuals numeric vector of OLS residuals
#' @param symmetry_threshold numeric: |gamma3| threshold for symmetry (default 0.3)
#' @param kurtosis_threshold numeric: gamma4 threshold for PMM3 (default -0.7)
#' @param g2_threshold numeric: minimum g2 improvement to justify PMM2 (default 0.95)
#' @param verbose logical: print decision reasoning (default TRUE)
#'
#' @return A list with components:
#'   \item{method}{Character: "OLS", "PMM2", or "PMM3"}
#'   \item{gamma3}{Sample skewness}
#'   \item{gamma4}{Sample excess kurtosis}
#'   \item{gamma6}{Sample 6th cumulant coefficient}
#'   \item{g2}{PMM2 efficiency factor}
#'   \item{g3}{PMM3 efficiency factor}
#'   \item{g_selected}{Efficiency factor for chosen method}
#'   \item{improvement_pct}{Expected variance reduction percentage}
#'   \item{reasoning}{Human-readable decision explanation}
#'   \item{n}{Sample size}
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(200); eps <- runif(200, -1, 1)
#' y <- 1 + 2 * x + eps
#' fit_ols <- lm(y ~ x)
#' pmm_dispatch(residuals(fit_ols))
#'
#' @export
pmm_dispatch <- function(residuals,
                         symmetry_threshold = 0.3,
                         kurtosis_threshold = -0.7,
                         g2_threshold       = 0.95,
                         verbose            = TRUE) {

  cum <- .residual_cumulants(residuals)
  g2  <- .g2_factor(cum$gamma3, cum$gamma4)
  g3  <- .g3_factor(cum$gamma4, cum$gamma6)

  abs_g3 <- abs(cum$gamma3)
  is_asymmetric          <- abs_g3 > symmetry_threshold
  is_strongly_asymmetric <- abs_g3 > 1.0
  is_platykurtic         <- cum$gamma4 < kurtosis_threshold

  # Decision logic
  if (is_strongly_asymmetric) {
    method <- "PMM2"
    g      <- g2
    reason <- sprintf(
      "|gamma3| = %.3f > 1.0: strong asymmetry detected. PMM2 expected to reduce variance by %.1f%%.",
      abs_g3, (1 - g2) * 100
    )
  } else if (is_asymmetric && g2 < g2_threshold) {
    method <- "PMM2"
    g      <- g2
    reason <- sprintf(
      "|gamma3| = %.3f > %.1f and g2 = %.3f < %.2f: moderate asymmetry, PMM2 worthwhile (%.1f%% reduction).",
      abs_g3, symmetry_threshold, g2, g2_threshold, (1 - g2) * 100
    )
  } else if (!is_asymmetric && is_platykurtic) {
    method <- "PMM3"
    g      <- g3
    reason <- sprintf(
      "|gamma3| = %.3f < %.1f (symmetric) and gamma4 = %.3f < %.1f (platykurtic). PMM3 expected to reduce variance by %.1f%%.",
      abs_g3, symmetry_threshold, cum$gamma4, kurtosis_threshold, (1 - g3) * 100
    )
  } else if (is_asymmetric && is_platykurtic) {
    method <- "PMM2"
    g      <- g2
    reason <- sprintf(
      "|gamma3| = %.3f (asymmetric) AND gamma4 = %.3f (platykurtic). Asymmetry dominates -- PMM2 (%.1f%% reduction).",
      abs_g3, cum$gamma4, (1 - g2) * 100
    )
  } else {
    method <- "OLS"
    g      <- 1.0
    reason <- sprintf(
      "gamma3 = %.3f, gamma4 = %.3f: near-Gaussian residuals. No PMM advantage expected. Use OLS.",
      cum$gamma3, cum$gamma4
    )
  }

  result <- list(
    method          = method,
    gamma3          = cum$gamma3,
    gamma4          = cum$gamma4,
    gamma6          = cum$gamma6,
    g2              = g2,
    g3              = g3,
    g_selected      = g,
    improvement_pct = (1 - g) * 100,
    reasoning       = reason,
    n               = length(residuals)
  )

  if (verbose) {
    cat(sprintf("  n = %d | gamma3 = %+.3f | gamma4 = %+.3f\n",
                result$n, result$gamma3, result$gamma4))
    cat(sprintf("  g2(PMM2) = %.4f | g3(PMM3) = %.4f\n", g2, g3))
    cat(sprintf("  >>> %s\n", reason))
  }

  invisible(result)
}
