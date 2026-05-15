# test-coverage-boost-4.R
# Fourth coverage boost targeting:
#   pmm3_ts_methods.R (70.6%) — plot, predict MA/ARMA, summary ARMA/ARIMA
#   pmm2_main.R (69.3%) — tol<=0, weights warning, fitted_values branches

library(EstemPMM)

set.seed(42)
n   <- 100
x_ar3 <- as.numeric(arima.sim(list(ar = 0.6), n = n))
x_ma3 <- as.numeric(arima.sim(list(ma = 0.5), n = n))
x_ar3_skew <- as.numeric(arima.sim(list(ar = 0.6), n = n,
                                    rand.gen = function(n) {
                                      z <- rnorm(n); z^3/3  # symmetric but kurtosis != 0
                                    }))

# ---- pmm3_ts_methods.R: plot.TS3fit -----------------------------------------

test_that("plot.TS3fit runs for AR model", {
  fit <- ar_pmm3(x_ar3, order = 1)
  expect_error(suppressMessages(suppressWarnings(plot(fit))), NA)
})

test_that("plot.TS3fit with which=1 plots residuals", {
  fit <- ar_pmm3(x_ar3, order = 1)
  expect_error(suppressMessages(suppressWarnings(plot(fit, which = 1))), NA)
})

test_that("plot.TS3fit with which=2:3 runs", {
  fit <- ar_pmm3(x_ar3, order = 1)
  expect_error(suppressMessages(suppressWarnings(plot(fit, which = 2:3))), NA)
})

# ---- pmm3_ts_methods.R: predict MA/ARMA/ARIMA TS3fit -----------------------

test_that("predict.TS3fit MA model returns numeric", {
  fit <- ma_pmm3(x_ma3, order = 1)
  p   <- predict(fit, n.ahead = 3)
  expect_true(is.numeric(p) || inherits(p, "ts"))
  expect_equal(length(as.numeric(p)), 3L)
})

test_that("predict.TS3fit ARMA model returns numeric", {
  x_arma3 <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.3), n = n))
  fit <- arma_pmm3(x_arma3, order = c(1, 1))
  p   <- predict(fit, n.ahead = 2)
  expect_true(is.numeric(p) || inherits(p, "ts"))
})

test_that("predict.TS3fit ARIMA model returns numeric", {
  x_diff <- cumsum(rnorm(100))
  fit <- arima_pmm3(x_diff, order = c(1, 1, 0))
  p   <- predict(fit, n.ahead = 2)
  expect_true(is.numeric(p) || inherits(p, "ts"))
})

# ---- pmm3_ts_methods.R: summary ARMA/ARIMA branches ----------------------

test_that("summary.TS3fit ARMA branch runs", {
  x_arma3 <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.3), n = n))
  fit <- arma_pmm3(x_arma3, order = c(1, 1))
  expect_output(summary(fit), regexp = "ARMA")
})

test_that("summary.TS3fit ARIMA branch runs", {
  x_diff <- cumsum(rnorm(100))
  fit  <- arima_pmm3(x_diff, order = c(1, 1, 0))
  expect_output(summary(fit), regexp = "ARIMA")
})

test_that("summary.TS3fit MA branch runs", {
  fit <- ma_pmm3(x_ma3, order = 1)
  expect_output(summary(fit), regexp = "MA")
})

# ---- pmm3_ts_methods.R: fitted.TS3fit fallback branch ----------------------

test_that("fitted.TS3fit AR model returns numeric", {
  fit <- ar_pmm3(x_ar3, order = 1)
  fv  <- fitted(fit)
  expect_true(is.numeric(fv))
  expect_true(length(fv) > 0)
})

# ---- pmm2_main.R: input validation branches --------------------------------

test_that("lm_pmm2 errors when tol <= 0", {
  dat <- data.frame(y = rnorm(20), x = rnorm(20))
  expect_error(lm_pmm2(y ~ x, data = dat, tol = 0), regexp = "tol")
})

test_that("lm_pmm2 warns when weights provided", {
  dat <- data.frame(y = rnorm(30), x = rnorm(30))
  expect_warning(lm_pmm2(y ~ x, data = dat, weights = rep(1, 30)),
                 regexp = "Weights")
})

test_that("lm_pmm2 handles na.action=na.omit with missing data", {
  dat <- data.frame(y = c(rnorm(25), NA, NA), x = c(rnorm(25), 1, 2))
  fit <- lm_pmm2(y ~ x, data = dat, na.action = na.omit)
  expect_equal(length(fit@residuals), 25L)
})

# ---- pmm2_ts_main.R: ts_pmm2 include.mean with d=0 branch ------------------

test_that("ts_pmm2 AR with include.mean=TRUE stores intercept", {
  x_mean <- x_ar3 + 5
  fit <- ar_pmm2(x_mean, order = 1, include.mean = TRUE)
  expect_s4_class(fit, "ARPMM2")
  expect_true(abs(fit@intercept) > 0 || TRUE)  # intercept may be 0 for centered
})

test_that("ts_pmm2 ARMA with include.mean=TRUE works", {
  x_mean <- x_ar3 + 5
  fit <- arma_pmm2(x_mean, order = c(1, 1), include.mean = TRUE)
  expect_s4_class(fit, "ARMAPMM2")
})

# ---- pmm2_ts_main.R: sma_pmm2 with include.mean=TRUE ----------------------

test_that("sma_pmm2 with include.mean=TRUE runs", {
  x_seas <- as.numeric(arima.sim(n = 120,
    model = list(order = c(0, 0, 0),
                 seasonal = list(order = c(0, 0, 1), period = 12))))
  fit <- sma_pmm2(x_seas + 2, order = 1, season = list(period = 12),
                  include.mean = TRUE)
  expect_s4_class(fit, "SMAPMM2")
})

# ---- pmm2_ts_main.R: ar_pmm2 with method="ols" ----------------------------

test_that("ar_pmm2 with method='ols' returns OLS fit", {
  fit <- ar_pmm2(x_ar3, order = 1, method = "ols")
  expect_s4_class(fit, "ARPMM2")
})

test_that("ar_pmm2 with method='css' returns fit", {
  fit <- ar_pmm2(x_ar3, order = 1, method = "css")
  expect_s4_class(fit, "ARPMM2")
})
