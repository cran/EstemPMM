set.seed(42)
n <- 150

test_that("ar_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ar = c(0.7, 0.2)), n = n)

  fit <- ar_pmm2(x, order = 2)

  expect_s4_class(fit, "ARPMM2")
  expect_length(fit@coefficients, 2)
  expect_gt(length(fit@residuals), 0)
})

test_that("ma_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ma = c(0.5, -0.3)), n = n)

  fit <- ma_pmm2(x, order = 2)

  expect_s4_class(fit, "MAPMM2")
  expect_length(fit@coefficients, 2)
})

test_that("arma_pmm2 produces valid fit object", {
  x <- arima.sim(model = list(ar = 0.7, ma = 0.3), n = n)

  fit <- arma_pmm2(x, order = c(1, 1))

  expect_s4_class(fit, "ARMAPMM2")
  expect_length(fit@coefficients, 2)
})

test_that("arima_pmm2 handles differencing", {
  x <- cumsum(rnorm(n))

  fit <- arima_pmm2(x, order = c(1, 1, 0))

  expect_s4_class(fit, "ARIMAPMM2")
  expect_length(fit@coefficients, 1)
})

test_that("ts_pmm2 dispatches to correct classes", {
  x <- arima.sim(model = list(ar = 0.7), n = n)

  fit_ar <- ts_pmm2(x, model_type = "ar", order = 1)
  expect_s4_class(fit_ar, "ARPMM2")

  fit_ma <- ts_pmm2(x, model_type = "ma", order = 1)
  expect_s4_class(fit_ma, "MAPMM2")
})

test_that("Time series residuals are numeric", {
  x <- arima.sim(model = list(ar = 0.7), n = n)
  fit <- ar_pmm2(x, order = 1)

  expect_true(is.numeric(fit@residuals))
  expect_false(any(is.na(fit@residuals)))
})

test_that("ar_pmm2 coefficients are bounded", {
  x <- arima.sim(model = list(ar = c(0.7, 0.2)), n = n)
  fit <- ar_pmm2(x, order = 2)

  expect_true(all(abs(fit@coefficients) < 2))
})
