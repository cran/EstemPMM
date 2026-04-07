# pmm3_utils.R - Utility functions for PMM3 moment estimation

#' Compute central moments for PMM3 from residuals
#'
#' Computes the second, fourth, and sixth central moments, along with
#' standardised cumulant coefficients gamma3, gamma4, gamma6, the theoretical
#' efficiency factor g3, and the moment ratio kappa used in the PMM3 solver.
#'
#' @param residuals numeric vector of residuals (typically from OLS)
#'
#' @return A list with components:
#'   \item{m2}{Second central moment}
#'   \item{m4}{Fourth central moment}
#'   \item{m6}{Sixth central moment}
#'   \item{gamma3}{Skewness coefficient (for symmetry check)}
#'   \item{gamma4}{Excess kurtosis}
#'   \item{gamma6}{Sixth-order cumulant coefficient}
#'   \item{g3}{Theoretical variance reduction factor}
#'   \item{kappa}{Moment ratio for NR solver (NA if near-Gaussian)}
#'   \item{improvement_pct}{Expected variance reduction percentage}
#'
#' @export
compute_moments_pmm3 <- function(residuals) {
  r  <- residuals - mean(residuals)
  m2 <- mean(r^2)
  m4 <- mean(r^4)
  m6 <- mean(r^6)

  gamma3 <- mean(r^3) / m2^(3/2)
  gamma4 <- m4 / m2^2 - 3
  gamma6 <- m6 / m2^3 - 15 * (m4 / m2^2) + 30

  denom_g3 <- 6 + 9 * gamma4 + gamma6
  g3 <- if (!is.na(denom_g3) && denom_g3 > 0) 1 - gamma4^2 / denom_g3 else 1.0

  denom_kappa <- m4 - 3 * m2^2
  kappa <- if (abs(denom_kappa) > .Machine$double.eps * 1e6) {
    (m6 - 3 * m4 * m2) / denom_kappa
  } else {
    NA_real_
  }

  list(m2 = m2, m4 = m4, m6 = m6,
       gamma3 = gamma3, gamma4 = gamma4, gamma6 = gamma6,
       g3 = g3, kappa = kappa,
       improvement_pct = (1 - g3) * 100)
}

#' Calculate PMM3 theoretical variance reduction factor
#'
#' Computes the standardised cumulant coefficients gamma4 and gamma6 from
#' the raw central moments, and derives the PMM3 efficiency factor
#' \eqn{g_3 = 1 - \gamma_4^2 / (6 + 9\gamma_4 + \gamma_6)}.
#'
#' @param m2 Second central moment
#' @param m4 Fourth central moment
#' @param m6 Sixth central moment
#'
#' @return A list with components:
#'   \item{gamma4}{Excess kurtosis}
#'   \item{gamma6}{Sixth-order cumulant coefficient}
#'   \item{g3}{Theoretical variance reduction factor}
#'
#' @export
pmm3_variance_factor <- function(m2, m4, m6) {
  if (is.na(m2) || m2 <= 0) {
    return(list(gamma4 = NA_real_, gamma6 = NA_real_, g3 = NA_real_))
  }
  gamma4 <- m4 / m2^2 - 3
  gamma6 <- m6 / m2^3 - 15 * (m4 / m2^2) + 30
  denom <- 6 + 9 * gamma4 + gamma6
  g3 <- if (is.na(denom) || denom <= 0) NA_real_ else 1 - gamma4^2 / denom
  list(gamma4 = gamma4, gamma6 = gamma6, g3 = g3)
}

#' Compute sixth-order cumulant coefficient gamma6
#'
#' Calculates \eqn{\gamma_6 = m_6/m_2^3 - 15 m_4/m_2^2 + 30} from a numeric
#' vector. For a Gaussian distribution gamma6 equals zero.
#'
#' @param x numeric vector
#'
#' @return Numeric scalar: the sixth-order cumulant coefficient
#' @export
pmm_gamma6 <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 6) {
    warning("At least 6 non-missing values are needed to compute gamma6")
    return(NA)
  }

  x_centered <- x - mean(x)
  m2 <- mean(x_centered^2)
  m4 <- mean(x_centered^4)
  m6 <- mean(x_centered^6)

  m6 / m2^3 - 15 * (m4 / m2^2) + 30
}

#' Test whether residuals are sufficiently symmetric for PMM3
#'
#' Computes the skewness coefficient gamma3 and checks whether its absolute
#' value falls below a given threshold. This helps decide between PMM2
#' (asymmetric) and PMM3 (symmetric platykurtic) estimation.
#'
#' @param x numeric vector of residuals
#' @param threshold numeric threshold for |gamma3| (default 0.3)
#'
#' @return A list with components:
#'   \item{gamma3}{Sample skewness coefficient}
#'   \item{is_symmetric}{Logical: TRUE if |gamma3| <= threshold}
#'   \item{message}{Human-readable verdict}
#'
#' @export
test_symmetry <- function(x, threshold = 0.3) {
  x <- x[!is.na(x)]
  x_centered <- x - mean(x)
  m2 <- mean(x_centered^2)
  m3 <- mean(x_centered^3)
  gamma3 <- m3 / m2^(3/2)
  is_sym <- abs(gamma3) <= threshold

  msg <- if (is_sym) {
    sprintf("|gamma3| = %.3f <= %.1f: residuals are symmetric. PMM3 may be appropriate.",
            abs(gamma3), threshold)
  } else {
    sprintf("|gamma3| = %.3f > %.1f: residuals are asymmetric. Consider PMM2 instead.",
            abs(gamma3), threshold)
  }

  list(gamma3 = gamma3, is_symmetric = is_sym, message = msg)
}
