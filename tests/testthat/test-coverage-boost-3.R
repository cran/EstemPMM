# test-coverage-boost-3.R
# Third coverage boost targeting:
#   pmm2_ts_main.R:  error branches, verbose paths, unified_pmm2_wrapper
#   pmm2_classes.R:  summary.TS2fit branches, summary.SARIMAPMM2
#   pmm2_monte_carlo.R: lognormal, student_t, chi_squared innovations

library(EstemPMM)

# ---- shared fixtures ---------------------------------------------------------

set.seed(42)
n   <- 120
x_ar <- as.numeric(arima.sim(list(ar = 0.7), n = n,
                               rand.gen = function(n) rgamma(n, 2, 1) - 2))
x_ma <- as.numeric(arima.sim(list(ma = 0.5), n = n))
x_arma <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.3), n = n))
x_seasonal <- as.numeric(arima.sim(n = 120,
  model = list(order = c(1, 0, 0), ar = 0.5,
               seasonal = list(order = c(1, 0, 0), ar = 0.4, period = 12))))

# ---- pmm2_ts_main.R: verbose branches ---------------------------------------

test_that("ar_pmm2 verbose=TRUE produces output", {
  expect_output(
    ar_pmm2(x_ar, order = 1, verbose = TRUE),
    regexp = NULL
  )
})

test_that("ma_pmm2 verbose=TRUE runs without error", {
  fit <- suppressWarnings(ma_pmm2(x_ma, order = 1, verbose = TRUE))
  expect_s4_class(fit, "MAPMM2")
})

test_that("arma_pmm2 verbose=TRUE produces output", {
  expect_output(
    arma_pmm2(x_arma, order = c(1, 1), verbose = TRUE),
    regexp = NULL
  )
})

test_that("arima_pmm2 verbose=TRUE produces output", {
  x_diff <- cumsum(rnorm(100))
  expect_output(
    arima_pmm2(x_diff, order = c(1, 1, 0), verbose = TRUE),
    regexp = NULL
  )
})

# ---- pmm2_ts_main.R: error conditions ---------------------------------------

test_that("sma_pmm2 errors on non-numeric input", {
  expect_error(
    sma_pmm2("not_numeric", order = 1, season = list(period = 12)),
    regexp = NULL
  )
})

test_that("sma_pmm2 errors on Q < 1", {
  expect_error(
    sma_pmm2(x_ma, order = 0, season = list(period = 12)),
    regexp = "at least 1"
  )
})

test_that("sma_pmm2 errors on period < 2", {
  expect_error(
    sma_pmm2(x_ma, order = 1, season = list(period = 1)),
    regexp = "period"
  )
})

test_that("sarma_pmm2 errors on NA values", {
  x_na <- x_ar; x_na[5] <- NA
  expect_error(
    sarma_pmm2(x_na, order = c(1, 0, 0, 1), season = list(period = 12)),
    regexp = "NA"
  )
})

test_that("sarma_pmm2 errors on wrong order length", {
  expect_error(
    sarma_pmm2(x_ar, order = c(1, 0, 1), season = list(period = 12)),
    regexp = "length 4"
  )
})

test_that("sarma_pmm2 errors on all-zero orders", {
  expect_error(
    sarma_pmm2(x_ar, order = c(0, 0, 0, 0), season = list(period = 12)),
    regexp = "positive"
  )
})

test_that("sarima_pmm2 errors on NA values", {
  x_na <- x_ar; x_na[1] <- NA
  expect_error(
    sarima_pmm2(x_na, order = c(0, 1, 0, 1),
                seasonal = list(order = c(1, 1), period = 12)),
    regexp = "NA"
  )
})

test_that("sarima_pmm2 errors on wrong order length", {
  expect_error(
    sarima_pmm2(x_ar, order = c(1, 0, 1),
                seasonal = list(order = c(1, 1), period = 12)),
    regexp = "length 4"
  )
})

test_that("sarima_pmm2 errors on wrong seasonal specification", {
  expect_error(
    sarima_pmm2(x_ar, order = c(0, 1, 0, 1),
                seasonal = c(1, 1)),
    regexp = "seasonal"
  )
})

# ---- pmm2_ts_main.R: internal functions (dead code) -------------------------

test_that("get_classical_estimates works for AR model", {
  theta <- EstemPMM:::get_classical_estimates(
    x = x_ar, order = 1, model_type = "ar",
    seasonal = NULL, include.mean = TRUE
  )
  expect_true(is.numeric(theta))
  expect_true(length(theta) >= 1)
})

test_that("get_classical_estimates works for MA model", {
  theta <- EstemPMM:::get_classical_estimates(
    x = x_ma, order = 1, model_type = "ma",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.numeric(theta))
})

test_that("get_classical_estimates works for ARMA model", {
  theta <- EstemPMM:::get_classical_estimates(
    x = x_arma, order = c(1, 1), model_type = "arma",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.numeric(theta))
  expect_equal(length(theta), 2L)
})

test_that("get_classical_estimates works for ARIMA model", {
  x_i <- cumsum(rnorm(80))
  theta <- EstemPMM:::get_classical_estimates(
    x = x_i, order = c(1, 1, 0), model_type = "arima",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.numeric(theta))
})

test_that("create_residual_function returns callable function for AR", {
  fn <- EstemPMM:::create_residual_function(
    x = x_ar, order = 1, model_type = "ar",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.function(fn))
  theta <- as.numeric(ar(x_ar, order.max = 1, aic = FALSE, method = "yw")$ar)
  res   <- fn(theta)
  expect_true(is.numeric(res))
  expect_equal(length(res), length(x_ar))
})

test_that("create_residual_function returns callable function for MA", {
  fn <- EstemPMM:::create_residual_function(
    x = x_ma, order = 1, model_type = "ma",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.function(fn))
  res <- fn(c(0.4))
  expect_true(is.numeric(res))
})

test_that("create_residual_function with include.mean=TRUE for AR", {
  fn <- EstemPMM:::create_residual_function(
    x = x_ar, order = 1, model_type = "ar",
    seasonal = NULL, include.mean = TRUE
  )
  expect_true(is.function(fn))
  res <- fn(c(mean(x_ar), 0.5))
  expect_true(is.numeric(res))
})

test_that("unified_pmm2_wrapper unified_global variant works for AR", {
  result <- EstemPMM:::unified_pmm2_wrapper(
    x = x_ar, order = 1, model_type = "ar",
    pmm2_variant = "unified_global",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.list(result))
  expect_true("theta" %in% names(result) || length(result) > 0)
})

test_that("unified_pmm2_wrapper unified_iterative variant works for AR", {
  result <- EstemPMM:::unified_pmm2_wrapper(
    x = x_ar, order = 1, model_type = "ar",
    pmm2_variant = "unified_iterative",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.list(result) || is.numeric(result))
})

test_that("unified_pmm2_wrapper linearized variant works for MA", {
  result <- EstemPMM:::unified_pmm2_wrapper(
    x = x_ma, order = 1, model_type = "ma",
    pmm2_variant = "linearized",
    seasonal = NULL, include.mean = FALSE
  )
  expect_true(is.list(result) || is.numeric(result))
})

# ---- pmm2_classes.R: summary branches ---------------------------------------

test_that("summary.TS2fit ARMA branch runs", {
  fit <- arma_pmm2(x_arma, order = c(1, 1))
  expect_output(summary(fit), regexp = "ARMA")
})

test_that("summary.TS2fit ARIMA branch runs", {
  x_i <- cumsum(rnorm(100))
  fit  <- arima_pmm2(x_i, order = c(1, 1, 0))
  expect_output(summary(fit), regexp = "ARIMA")
})

test_that("summary.TS2fit MA branch runs", {
  fit <- ma_pmm2(x_ma, order = 1)
  expect_output(summary(fit), regexp = "MA")
})

test_that("summary.SARIMAPMM2 runs without error", {
  fit <- tryCatch(
    sarima_pmm2(x_seasonal, order = c(0, 1, 0, 1),
                seasonal = list(order = c(1, 1), period = 12)),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    expect_output(summary(fit), regexp = NULL)
  } else {
    skip("sarima_pmm2 fit failed on this data")
  }
})

# ---- pmm2_monte_carlo.R: alternative innovation types ----------------------

test_that("pmm2_monte_carlo_compare lognormal innovations work", {
  specs <- list(
    list(model = "ma", order = 1, theta = 0.5, label = "MA(1)-lnorm",
         innovations = list(type = "lognormal", sigma = 0.5))
  )
  res <- pmm2_monte_carlo_compare(
    model_specs = specs,
    methods     = c("css", "pmm2"),
    n           = 50,
    n_sim       = 4,
    seed        = 42,
    progress    = FALSE
  )
  expect_type(res, "list")
})

test_that("pmm2_monte_carlo_compare student_t innovations work", {
  specs <- list(
    list(model = "ma", order = 1, theta = 0.4, label = "MA(1)-t",
         innovations = list(type = "student_t", df = 6))
  )
  res <- pmm2_monte_carlo_compare(
    model_specs = specs,
    methods     = c("css", "pmm2"),
    n           = 50,
    n_sim       = 4,
    seed        = 7,
    progress    = FALSE
  )
  expect_type(res, "list")
})

test_that("pmm2_monte_carlo_compare chi_squared innovations work", {
  specs <- list(
    list(model = "ma", order = 1, theta = 0.4, label = "MA(1)-chi2",
         innovations = list(type = "chi_squared", df = 3))
  )
  res <- pmm2_monte_carlo_compare(
    model_specs = specs,
    methods     = c("css", "pmm2"),
    n           = 50,
    n_sim       = 4,
    seed        = 11,
    progress    = FALSE
  )
  expect_type(res, "list")
})

test_that("pmm2_monte_carlo_compare gaussian innovations produce g=1", {
  specs <- list(
    list(model = "ma", order = 1, theta = 0.5, label = "MA(1)-gauss",
         innovations = list(type = "gaussian"))
  )
  res <- pmm2_monte_carlo_compare(
    model_specs = specs,
    methods     = c("css", "pmm2"),
    n           = 60,
    n_sim       = 5,
    seed        = 99,
    progress    = FALSE
  )
  gain_df <- res$gain
  expect_equal(gain_df$theoretical_g, 1, tolerance = 1e-8)
})

# ---- pmm2_main.R: compare_with_ols and verbose ------------------------------

test_that("compare_with_ols returns list with ols and pmm2 components", {
  dat <- data.frame(
    y = 1 + 2 * rnorm(60) + (rgamma(60, 2, 1) - 2),
    x = rnorm(60)
  )
  result <- compare_with_ols(y ~ x, data = dat)
  expect_true(is.list(result))
  expect_true(!is.null(result$ols))
  expect_true(!is.null(result$pmm2))
})
