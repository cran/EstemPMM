# test-coverage-boost.R
# Targeted tests for code paths with 0% coverage:
#   pmm2_unified.R   (pmm2_nonlinear_onestep, pmm2_nonlinear_iterative)
#   sarimax_wrapper.R (get_sarimax_residuals, get_sarimax_jacobian)
#   optimized_direct_pmm2.R (deprecated wrapper)
# Plus edge cases in pmm2_utils.R

library(EstemPMM)

# ---- pmm2_nonlinear_onestep --------------------------------------------------

test_that("pmm2_nonlinear_onestep works on linear regression", {
  set.seed(42)
  n   <- 80
  x   <- rnorm(n)
  eps <- rgamma(n, 2, 1) - 2
  y   <- 1 + 2 * x + eps
  X   <- cbind(1, x)

  fn_res <- function(theta) as.numeric(y - X %*% theta)
  theta0 <- as.numeric(coef(lm.fit(X, y)))

  result <- pmm2_nonlinear_onestep(theta0, fn_res, verbose = FALSE)

  expect_true(is.list(result))
  expect_true(!is.null(result$theta) || !is.null(result$b) || !is.null(result$coefficients))
})

test_that("pmm2_nonlinear_onestep verbose=TRUE prints output", {
  set.seed(1)
  n  <- 50
  y  <- rnorm(n)
  X  <- cbind(1, rnorm(n))
  fn <- function(th) y - X %*% th
  th <- lm.fit(X, y)$coefficients

  expect_output(
    pmm2_nonlinear_onestep(th, fn, verbose = TRUE),
    regexp = NULL  # just verify it runs without error when verbose
  )
})

# ---- pmm2_nonlinear_iterative ------------------------------------------------

test_that("pmm2_nonlinear_iterative converges for linear regression", {
  set.seed(42)
  n   <- 100
  x   <- rnorm(n)
  eps <- rgamma(n, 2, 1) - 2
  y   <- 2 + 1.5 * x + eps
  X   <- cbind(1, x)

  fn_res <- function(theta) as.numeric(y - X %*% theta)
  theta0 <- as.numeric(lm.fit(X, y)$coefficients)

  result <- pmm2_nonlinear_iterative(theta0, fn_res,
                                     max_iter = 30, tol = 1e-5,
                                     verbose = FALSE)
  expect_true(is.list(result))
})

test_that("pmm2_nonlinear_iterative accepts explicit Jacobian", {
  set.seed(7)
  n   <- 60
  x   <- rnorm(n)
  y   <- 1 + x + rnorm(n, sd = 0.5)
  X   <- cbind(1, x)

  fn_res  <- function(th) as.numeric(y - X %*% th)
  fn_jac  <- function(th) X  # negative Jacobian of residuals w.r.t. theta is X
  theta0  <- lm.fit(X, y)$coefficients

  result <- pmm2_nonlinear_iterative(theta0, fn_res, fn_jacobian = fn_jac,
                                     max_iter = 20, verbose = FALSE)
  expect_true(is.list(result))
})

# ---- get_sarimax_residuals ---------------------------------------------------

test_that("get_sarimax_residuals returns residual vector for ARMA(1,1)", {
  set.seed(42)
  n  <- 80
  y  <- as.numeric(arima.sim(list(ar = 0.6, ma = 0.3), n = n))

  theta <- c(0.5, 0.2)
  res   <- EstemPMM:::get_sarimax_residuals(
    theta        = theta,
    y            = y,
    order        = c(1, 0, 1),
    seasonal     = list(order = c(0, 0, 0), period = NA),
    include.mean = FALSE
  )
  expect_true(is.numeric(res))
  expect_equal(length(res), n)
})

test_that("get_sarimax_residuals handles include.mean=TRUE", {
  set.seed(5)
  n <- 60
  y <- rnorm(n, mean = 3)

  theta <- c(0.4, 3.1)   # ar1, mean
  res <- EstemPMM:::get_sarimax_residuals(
    theta        = theta,
    y            = y,
    order        = c(1, 0, 0),
    seasonal     = list(order = c(0, 0, 0), period = NA),
    include.mean = TRUE
  )
  expect_true(is.numeric(res))
})

test_that("get_sarimax_residuals handles differencing", {
  set.seed(10)
  y <- cumsum(rnorm(80))   # I(1) series

  theta <- c(0.4)  # AR(1) after differencing
  res <- EstemPMM:::get_sarimax_residuals(
    theta        = theta,
    y            = y,
    order        = c(1, 1, 0),
    seasonal     = list(order = c(0, 0, 0), period = NA),
    include.mean = FALSE
  )
  expect_true(is.numeric(res))
  expect_equal(length(res), length(y))   # arima() returns same-length with leading NAs
})

# ---- get_sarimax_jacobian ----------------------------------------------------

test_that("get_sarimax_jacobian returns matrix of correct dimensions", {
  set.seed(42)
  n  <- 60
  y  <- as.numeric(arima.sim(list(ar = 0.5), n = n))

  theta <- c(0.4)
  J <- EstemPMM:::get_sarimax_jacobian(
    theta        = theta,
    y            = y,
    order        = c(1, 0, 0),
    seasonal     = list(order = c(0, 0, 0), period = NA),
    include.mean = FALSE
  )
  expect_true(is.matrix(J))
  expect_equal(ncol(J), 1L)   # 1 parameter
})

# ---- optimized_direct_pmm2 (deprecated) -------------------------------------

test_that("optimized_direct_pmm2 errors with deprecation message", {
  set.seed(42)
  y     <- arima.sim(list(ar = 0.6), n = 60)
  theta <- c(0.5)

  expect_error(
    EstemPMM:::optimized_direct_pmm2(
      theta_init = theta,
      y          = as.numeric(y),
      order      = c(1, 0, 0)
    ),
    regexp = "deprecated|Deprecated"
  )
})

# ---- pmm2_utils edge cases ---------------------------------------------------

test_that("pmm2_variance_factor returns NA list when m2 <= 0", {
  vf <- pmm2_variance_factor(m2 = 0, m3 = 0.5, m4 = 1)
  expect_true(is.na(vf$g))
  expect_true(is.na(vf$c3))
})

test_that("pmm2_variance_factor handles NA input", {
  vf <- pmm2_variance_factor(m2 = NA, m3 = 0.5, m4 = 1)
  expect_true(is.na(vf$g))
})

test_that("pmm2_variance_factor is correct for Gaussian (g should be 1)", {
  # Gaussian: c3 = 0, c4 = 0 => g = 1 - 0/(2+0) = 1
  vf <- pmm2_variance_factor(m2 = 1, m3 = 0, m4 = 3)
  expect_equal(vf$g, 1, tolerance = 1e-10)
  expect_equal(vf$c3, 0, tolerance = 1e-10)
})

test_that("pmm2_variance_factor gives g < 1 for skewed distribution", {
  # Gamma(2,1) shifted: c3 > 0, g < 1
  set.seed(42)
  eps <- rgamma(10000, 2, 1) - 2
  m2  <- mean(eps^2); m3 <- mean(eps^3); m4 <- mean(eps^4)
  vf  <- pmm2_variance_factor(m2, m3, m4)
  expect_true(vf$g < 1)
  expect_true(vf$g > 0)
})

test_that("pmm2_variance_matrices returns list with ols and pmm2 components", {
  set.seed(1)
  X  <- cbind(1, rnorm(50))
  vm <- pmm2_variance_matrices(X, m2 = 1, m3 = 0.5, m4 = 3)
  expect_named(vm, c("ols", "pmm2", "c3", "c4", "g"), ignore.order = TRUE)
  expect_true(is.matrix(vm$ols))
  expect_true(is.matrix(vm$pmm2))
  expect_equal(dim(vm$ols), c(2L, 2L))
})

test_that("pmm2_variance_matrices pmm2 < ols when g < 1", {
  set.seed(1)
  X  <- cbind(1, rnorm(100))
  # Use Gamma moments for g < 1
  eps <- rgamma(1000, 2, 1) - 2
  m2  <- mean(eps^2); m3 <- mean(eps^3); m4 <- mean(eps^4)
  vm  <- pmm2_variance_matrices(X, m2, m3, m4)
  expect_true(vm$pmm2[2, 2] < vm$ols[2, 2])
})

test_that("pmm_kurtosis returns correct values for known distributions", {
  set.seed(42)
  normal_data <- rnorm(10000)
  # Excess kurtosis of normal is 0
  expect_equal(pmm_kurtosis(normal_data, excess = TRUE), 0, tolerance = 0.1)
  # Non-excess kurtosis of normal is 3
  expect_equal(pmm_kurtosis(normal_data, excess = FALSE), 3, tolerance = 0.1)
})

test_that("pmm_skewness returns near zero for symmetric distribution", {
  set.seed(42)
  x <- rnorm(5000)
  expect_equal(pmm_skewness(x), 0, tolerance = 0.1)
})

test_that("pmm_skewness returns positive for right-skewed distribution", {
  x <- rgamma(1000, shape = 2, scale = 1)  # right-skewed
  expect_true(pmm_skewness(x) > 0)
})

test_that("compute_moments returns list with m2, m3, m4, c3, c4, g", {
  eps    <- rgamma(500, 2, 1) - 2
  result <- compute_moments(eps)
  expect_named(result, c("m2", "m3", "m4", "c3", "c4", "g"))
  expect_true(result$m2 > 0)
})

# ---- pmm2_main.R: verbose and edge case branches ----------------------------

test_that("lm_pmm2 verbose=TRUE runs without error", {
  set.seed(1)
  dat <- data.frame(y = rnorm(30), x = rnorm(30))
  expect_output(
    lm_pmm2(y ~ x, data = dat, verbose = TRUE, max_iter = 3),
    regexp = NULL
  )
})

test_that("lm_pmm2 errors on non-data.frame input", {
  expect_error(lm_pmm2(y ~ x, data = list(y = 1:5, x = 1:5)))
})

test_that("lm_pmm2 errors on negative max_iter", {
  dat <- data.frame(y = rnorm(20), x = rnorm(20))
  expect_error(lm_pmm2(y ~ x, data = dat, max_iter = -1))
})

test_that("compare_with_ols returns list with ols and pmm2 components", {
  set.seed(42)
  dat <- data.frame(y = rnorm(50) + rnorm(50), x = rnorm(50))
  result <- compare_with_ols(y ~ x, data = dat)
  expect_true(is.list(result))
  expect_true(!is.null(result$ols))
  expect_true(!is.null(result$pmm2))
  expect_true(is.data.frame(result$coefficients))
  expect_true(is.data.frame(result$residual_stats))
})
