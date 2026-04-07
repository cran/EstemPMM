# Tests for PMM3 time series models

set.seed(42)
n <- 200

# Generate AR(1) with uniform innovations (platykurtic)
innov_unif <- runif(n + 50, -sqrt(3), sqrt(3))
ar_data <- arima.sim(n = n, list(ar = 0.7), innov = innov_unif[1:n])

# Generate MA(1) with uniform innovations
ma_data <- arima.sim(n = n, list(ma = 0.5), innov = innov_unif[1:n])

# =============================================================================
# AR models
# =============================================================================

test_that("ar_pmm3 produces valid fit", {
  fit <- ar_pmm3(ar_data, order = 1)

  expect_s4_class(fit, "ARPMM3")
  expect_true(fit@convergence)
  expect_length(fit@coefficients, 1)
  expect_true(fit@m2 > 0)
  expect_true(fit@m4 > 0)
  expect_true(fit@m6 > 0)
})

test_that("ar_pmm3 coefficients are reasonable", {
  fit <- ar_pmm3(ar_data, order = 1)

  expect_true(is.numeric(coef(fit)))
  expect_false(any(is.na(fit@coefficients)))
  # AR(1) coef should be near 0.7
  expect_true(abs(fit@coefficients[1] - 0.7) < 0.3)
})

test_that("ar_pmm3 S4 methods work", {
  fit <- ar_pmm3(ar_data, order = 1)

  expect_true(length(coef(fit)) >= 1)
  expect_true(length(residuals(fit)) > 0)
  expect_true(length(fitted(fit)) > 0)
  expect_output(summary(fit), "PMM3 time series")
})

test_that("ar_pmm3 with higher order works", {
  fit <- ar_pmm3(ar_data, order = 2)

  expect_s4_class(fit, "ARPMM3")
  expect_length(fit@coefficients, 2)
})

test_that("ar_pmm3 predict works", {
  fit <- ar_pmm3(ar_data, order = 1)
  preds <- predict(fit, n.ahead = 3)

  expect_length(preds, 3)
  expect_true(is.numeric(preds))
  expect_false(any(is.na(preds)))
})

test_that("ar_pmm3 adaptive mode works", {
  fit <- ar_pmm3(ar_data, order = 1, adaptive = TRUE)

  expect_s4_class(fit, "ARPMM3")
  expect_true(fit@convergence)
})

# =============================================================================
# MA models
# =============================================================================

test_that("ma_pmm3 produces valid fit", {
  fit <- ma_pmm3(ma_data, order = 1)

  expect_s4_class(fit, "MAPMM3")
  expect_true(fit@convergence)
  expect_length(fit@coefficients, 1)
})

test_that("ma_pmm3 S4 methods work", {
  fit <- ma_pmm3(ma_data, order = 1)

  expect_true(length(coef(fit)) >= 1)
  expect_true(length(residuals(fit)) > 0)
  expect_output(summary(fit), "PMM3 time series")
})

# =============================================================================
# ARMA models
# =============================================================================

test_that("arma_pmm3 produces valid fit", {
  arma_data <- arima.sim(n = n, list(ar = 0.5, ma = 0.3),
                         innov = innov_unif[1:n])
  fit <- arma_pmm3(arma_data, order = c(1, 1))

  expect_s4_class(fit, "ARMAPMM3")
  expect_true(fit@convergence)
  expect_length(fit@coefficients, 2)
})

# =============================================================================
# ARIMA models
# =============================================================================

test_that("arima_pmm3 produces valid fit", {
  arima_data <- cumsum(arima.sim(n = n, list(ar = 0.5),
                                 innov = innov_unif[1:n]))
  fit <- arima_pmm3(arima_data, order = c(1, 1, 0))

  expect_s4_class(fit, "ARIMAPMM3")
  expect_true(fit@convergence)
})

# =============================================================================
# General TS3fit checks
# =============================================================================

test_that("AIC method works for TS3fit", {
  fit <- ar_pmm3(ar_data, order = 1)
  aic <- AIC(fit)

  expect_true(is.numeric(aic))
  expect_true(is.finite(aic))
})

test_that("ts_pmm3 warns for asymmetric innovations", {
  set.seed(99)
  innov_exp <- rexp(n) - 1
  asym_data <- arima.sim(n = n, list(ar = 0.5), innov = innov_exp)

  expect_warning(ar_pmm3(asym_data, order = 1), "asymmetric|PMM2")
})

test_that("PMM3 moment slots are populated for TS3fit", {
  fit <- ar_pmm3(ar_data, order = 1)

  expect_true(!is.na(fit@gamma4))
  expect_true(!is.na(fit@gamma6))
  expect_true(!is.na(fit@g_coefficient))
  expect_true(fit@g_coefficient <= 1)
})
