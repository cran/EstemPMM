# Tests for Seasonal Models: SAR, SMA, SARMA, SARIMA
# Testing new functionality added in version 0.1.2+

set.seed(123)

# =============================================================================
# SAR (Seasonal AutoRegressive) Model Tests
# =============================================================================

test_that("sar_pmm2 produces valid SARPMM2 fit object", {
  # Simulate seasonal AR data: SAR(1)_12
  n <- 120
  phi_s <- 0.6
  x <- numeric(n)
  innovations <- rnorm(n, sd = 1)

  for (t in 13:n) {
    x[t] <- phi_s * x[t - 12] + innovations[t]
  }

  # Fit SAR model
  fit <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))

  # Check class
  expect_s4_class(fit, "SARPMM2")

  # Check coefficients
  expect_length(fit@coefficients, 1)  # Only 1 seasonal AR coefficient
  expect_true(is.numeric(fit@coefficients))
  expect_false(any(is.na(fit@coefficients)))

  # Check residuals
  expect_true(is.numeric(fit@residuals))
  expect_gt(length(fit@residuals), 0)

  # Check convergence
  expect_true(is.logical(fit@convergence))

  # Check model type
  expect_equal(fit@model_type, "sar")
})

test_that("sar_pmm2 handles multiplicative SAR models", {
  n <- 120
  phi <- 0.5  # Non-seasonal AR
  phi_s <- 0.4  # Seasonal AR
  x <- numeric(n)
  innovations <- rnorm(n, sd = 1)

  for (t in 14:n) {
    x[t] <- phi * x[t - 1] + phi_s * x[t - 12] - phi * phi_s * x[t - 13] + innovations[t]
  }

  # Fit multiplicative SAR(1,1)_12
  fit <- sar_pmm2(x, order = c(1, 1), season = list(period = 12),
                  multiplicative = TRUE)

  expect_s4_class(fit, "SARPMM2")
  expect_equal(length(fit@coefficients), 3)  # AR(1) + SAR(1) + Interaction
})

test_that("sar_pmm2 coefficients are bounded", {
  n <- 100
  x <- rnorm(n)

  fit <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))

  # Coefficients should be reasonable
  expect_true(all(abs(fit@coefficients) < 2))
})

test_that("sar_pmm2 works with different seasonal periods", {
  n <- 60
  x <- rnorm(n)

  # Quarterly data (period = 4)
  fit_q <- sar_pmm2(x, order = c(0, 1), season = list(period = 4))
  expect_s4_class(fit_q, "SARPMM2")

  # Monthly data (period = 12) - already tested above
  # Weekly data (period = 52) would need more data
})

test_that("sar_pmm2 handles mean/intercept correctly", {
  n <- 100
  x <- rnorm(n, mean = 5)

  fit_with_mean <- sar_pmm2(x, order = c(0, 1), season = list(period = 12),
                            include.mean = TRUE)
  fit_no_mean <- sar_pmm2(x, order = c(0, 1), season = list(period = 12),
                          include.mean = FALSE)

  expect_s4_class(fit_with_mean, "SARPMM2")
  expect_s4_class(fit_no_mean, "SARPMM2")

  # Intercept should differ
  expect_true(abs(fit_with_mean@intercept - fit_no_mean@intercept) > 0.1)
})

# =============================================================================
# SMA (Seasonal Moving Average) Model Tests
# =============================================================================

test_that("sma_pmm2 produces valid SMAPMM2 fit object", {
  n <- 120
  theta_s <- 0.5
  innovations <- rnorm(n, sd = 1)
  x <- numeric(n)

  for (t in 1:n) {
    if (t > 12) {
      x[t] <- innovations[t] + theta_s * innovations[t - 12]
    } else {
      x[t] <- innovations[t]
    }
  }

  # Fit SMA model
  fit <- sma_pmm2(x, order = 1, season = list(period = 12))

  # Check class
  expect_s4_class(fit, "SMAPMM2")

  # Check coefficients
  expect_length(fit@coefficients, 1)  # One seasonal MA coefficient
  expect_true(is.numeric(fit@coefficients))

  # Check innovations slot (specific to SMAPMM2)
  expect_true(is.numeric(fit@innovations))
  expect_equal(length(fit@innovations), n)

  # Check model type
  expect_equal(fit@model_type, "sma")
})

test_that("sma_pmm2 handles higher order models", {
  n <- 120
  x <- rnorm(n)

  # SMA(2)_12
  fit <- sma_pmm2(x, order = 2, season = list(period = 12))

  expect_s4_class(fit, "SMAPMM2")
  expect_length(fit@coefficients, 2)
})

test_that("sma_pmm2 handles different methods", {
  n <- 120
  x <- rnorm(n)

  # PMM2 method
  fit_pmm2 <- sma_pmm2(x, order = 1, season = list(period = 12), method = "pmm2")
  expect_s4_class(fit_pmm2, "SMAPMM2")

  # CSS method
  fit_css <- sma_pmm2(x, order = 1, season = list(period = 12), method = "css")
  expect_s4_class(fit_css, "SMAPMM2")

  # Coefficients should be similar but not identical
  expect_true(abs(fit_pmm2@coefficients[1] - fit_css@coefficients[1]) < 0.5)
})

test_that("sma_pmm2 convergence parameters work", {
  n <- 100
  x <- rnorm(n)

  fit <- sma_pmm2(x, order = 1, season = list(period = 12),
                  max_iter = 100, tol = 1e-6)

  expect_s4_class(fit, "SMAPMM2")
  expect_true(fit@iterations <= 100)
})

# =============================================================================
# SARMA (Seasonal ARMA) Model Tests
# =============================================================================

test_that("sarma_pmm2 produces valid SARMAPMM2 fit object", {
  n <- 120
  x <- rnorm(n)

  # SARMA(0,1,0,1)_12 - pure seasonal model
  fit <- sarma_pmm2(x, order = c(0, 1, 0, 1), season = list(period = 12))

  expect_s4_class(fit, "SARMAPMM2")
  expect_length(fit@coefficients, 2)  # 1 SAR + 1 SMA
  expect_equal(fit@model_type, "sarma")
})

test_that("sarma_pmm2 handles full SARMA(p,P,q,Q) models", {
  n <- 150
  x <- rnorm(n)

  # SARMA(1,1,1,1)_12
  fit <- sarma_pmm2(x, order = c(1, 1, 1, 1), season = list(period = 12))

  expect_s4_class(fit, "SARMAPMM2")
  # Should have: 1 AR + 1 SAR + 1 MA + 1 SMA = 4 coefficients
  expect_equal(length(fit@coefficients), 4)
})

test_that("sarma_pmm2 order slot contains correct information", {
  n <- 120
  x <- rnorm(n)

  fit <- sarma_pmm2(x, order = c(1, 1, 1, 1), season = list(period = 12))

  expect_true(is.list(fit@order))
  expect_true("period" %in% names(fit@order))
  expect_equal(fit@order$period, 12)
})

test_that("sarma_pmm2 handles mean parameter", {
  n <- 120
  x <- rnorm(n, mean = 3)

  fit_mean <- sarma_pmm2(x, order = c(0, 1, 0, 1), season = list(period = 12),
                         include.mean = TRUE)

  expect_s4_class(fit_mean, "SARMAPMM2")
  expect_true(abs(fit_mean@intercept) > 0.1)
})

test_that("sarma_pmm2 residuals have correct properties", {
  n <- 120
  x <- rnorm(n)

  fit <- sarma_pmm2(x, order = c(0, 1, 0, 1), season = list(period = 12))

  expect_true(is.numeric(fit@residuals))
  expect_false(any(is.na(fit@residuals)))
  expect_equal(length(fit@residuals), n)

  # Mean of residuals should be close to 0
  expect_true(abs(mean(fit@residuals)) < 0.5)
})

# =============================================================================
# SARIMA (Seasonal ARIMA with differencing) Model Tests
# =============================================================================

test_that("sarima_pmm2 produces valid SARIMAPMM2 fit object", {
  n <- 120
  # Create non-stationary seasonal data
  x <- cumsum(rnorm(n))

  # SARIMA(0,1,0,1) x (0,1,1,1)_12
  fit <- sarima_pmm2(x, order = c(0, 1, 0, 1),
                     seasonal = list(order = c(0, 1), period = 12))

  expect_s4_class(fit, "SARIMAPMM2")
  expect_equal(fit@model_type, "sarima")
})

test_that("sarima_pmm2 handles differencing orders", {
  n <- 120
  x <- cumsum(cumsum(rnorm(n)))  # Needs d=2

  fit <- sarima_pmm2(x, order = c(0, 2, 0, 1),
                     seasonal = list(order = c(0, 1), period = 12))

  expect_s4_class(fit, "SARIMAPMM2")
  expect_true("d" %in% names(fit@order))
  expect_true("D" %in% names(fit@order))
})

test_that("sarima_pmm2 full model specification", {
  n <- 150
  x <- rnorm(n)

  # SARIMA(1,1,1,1) x (1,1,1,1)_12
  fit <- sarima_pmm2(x, order = c(1, 1, 1, 1),
                     seasonal = list(order = c(1, 1), period = 12))

  expect_s4_class(fit, "SARIMAPMM2")
  # AR + SAR + MA + SMA coefficients
  expect_true(length(fit@coefficients) > 0)
})

test_that("sarima_pmm2 handles mean parameter appropriately", {
  n <- 120
  x <- rnorm(n)

  # With differencing, mean is typically NULL or auto-determined
  fit <- sarima_pmm2(x, order = c(0, 1, 0, 1),
                     seasonal = list(order = c(0, 1), period = 12),
                     include.mean = NULL)

  expect_s4_class(fit, "SARIMAPMM2")
})

# =============================================================================
# S4 Class Methods Tests
# =============================================================================

test_that("coef method works for SARPMM2", {
  n <- 100
  x <- rnorm(n)
  fit <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))

  coefs <- coef(fit)
  expect_true(is.numeric(coefs))
  expect_equal(length(coefs), length(fit@coefficients))
})

test_that("coef method works for SMAPMM2", {
  n <- 100
  x <- rnorm(n)
  fit <- sma_pmm2(x, order = 1, season = list(period = 12))

  coefs <- coef(fit)
  expect_true(is.numeric(coefs))
  expect_equal(length(coefs), 1)
})

test_that("residuals method works for seasonal models", {
  n <- 100
  x <- rnorm(n)

  fit_sar <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))
  resid_sar <- residuals(fit_sar)
  expect_true(is.numeric(resid_sar))
  expect_equal(length(resid_sar), n)

  fit_sma <- sma_pmm2(x, order = 1, season = list(period = 12))
  resid_sma <- residuals(fit_sma)
  expect_true(is.numeric(resid_sma))
})

test_that("fitted method works for seasonal models", {
  n <- 100
  x <- rnorm(n)

  fit_sar <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))
  fitted_vals <- fitted(fit_sar)

  expect_true(is.numeric(fitted_vals))
  # Fitted values + residuals should approximately equal original
  expect_true(max(abs(fitted_vals + residuals(fit_sar) - x)) < 0.1)
})

test_that("summary method works for SARPMM2", {
  n <- 100
  x <- rnorm(n)
  fit <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))

  # Summary should produce output
  expect_output(summary(fit))
})

test_that("summary method works for SMAPMM2", {
  n <- 100
  x <- rnorm(n)
  fit <- sma_pmm2(x, order = 1, season = list(period = 12))

  expect_output(summary(fit))
})

test_that("summary method works for SARMAPMM2", {
  n <- 120
  x <- rnorm(n)
  fit <- sarma_pmm2(x, order = c(0, 1, 0, 1), season = list(period = 12))

  expect_output(summary(fit))
})

test_that("summary method works for SARIMAPMM2", {
  n <- 120
  x <- cumsum(rnorm(n))
  fit <- sarima_pmm2(x, order = c(0, 1, 0, 1),
                     seasonal = list(order = c(0, 1), period = 12))

  expect_output(summary(fit))
})

# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

test_that("seasonal models handle short time series gracefully", {
  n <- 30  # Barely enough for period=12
  x <- rnorm(n)

  # Should work but may have convergence issues
  fit <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))
  expect_s4_class(fit, "SARPMM2")
})

test_that("seasonal models validate period parameter", {
  n <- 100
  x <- rnorm(n)

  # Very large period relative to data length might cause issues
  # but should not crash
  expect_s4_class(
    sar_pmm2(x, order = c(0, 1), season = list(period = 4)),
    "SARPMM2"
  )
})

test_that("seasonal models handle constant series", {
  n <- 100
  x <- rep(5, n)

  # Should handle constant series without crashing
  expect_error(
    sar_pmm2(x, order = c(0, 1), season = list(period = 12)),
    NA  # NA means we expect no error, or we expect it to work
  )
})

# =============================================================================
# Comparison with Base Stats Models
# =============================================================================

test_that("SAR model coefficients are reasonable compared to arima()", {
  n <- 120
  phi_s <- 0.6
  x <- numeric(n)
  innovations <- rnorm(n, sd = 1)

  for (t in 13:n) {
    x[t] <- phi_s * x[t - 12] + innovations[t]
  }

  fit_pmm2 <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))

  # Base R arima for comparison (if available)
  if (requireNamespace("stats", quietly = TRUE)) {
    fit_arima <- try(
      stats::arima(x, order = c(0, 0, 0),
                   seasonal = list(order = c(1, 0, 0), period = 12)),
      silent = TRUE
    )

    if (!inherits(fit_arima, "try-error")) {
      # Coefficients should be in the same ballpark
      expect_true(abs(fit_pmm2@coefficients[1] - coef(fit_arima)[1]) < 0.3)
    }
  }
})

# =============================================================================
# Integration Tests
# =============================================================================

test_that("All seasonal model classes inherit from TS2fit", {
  n <- 100
  x <- rnorm(n)

  fit_sar <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))
  fit_sma <- sma_pmm2(x, order = 1, season = list(period = 12))
  fit_sarma <- sarma_pmm2(x, order = c(0, 1, 0, 1), season = list(period = 12))

  expect_true(methods::is(fit_sar, "TS2fit"))
  expect_true(methods::is(fit_sma, "TS2fit"))
  expect_true(methods::is(fit_sarma, "TS2fit"))
})

test_that("Seasonal models store original series", {
  n <- 100
  x <- rnorm(n)

  fit <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))

  expect_equal(length(fit@original_series), n)
  expect_equal(fit@original_series, x)
})

test_that("Moment statistics are computed for seasonal models", {
  n <- 120
  x <- rnorm(n)

  fit <- sar_pmm2(x, order = c(0, 1), season = list(period = 12))

  expect_true(is.numeric(fit@m2))
  expect_true(is.numeric(fit@m3))
  expect_true(is.numeric(fit@m4))
  expect_true(fit@m2 > 0)  # Variance should be positive
})
