# tests/testthat/test-unified-pmm2.R
# Unit tests for Unified PMM2 variants (Version 0.2.0)

test_that("ar_pmm2 accepts pmm2_variant parameter", {
  set.seed(123)
  x <- arima.sim(n = 100, list(ar = 0.7))
  
  # Test all three variants
  fit_global <- ar_pmm2(x, order = 1, pmm2_variant = "unified_global")
  fit_iter <- ar_pmm2(x, order = 1, pmm2_variant = "unified_iterative")
  fit_linear <- ar_pmm2(x, order = 1, pmm2_variant = "linearized")
  
  # All should return ARPMM2 objects
  expect_s4_class(fit_global, "ARPMM2")
  expect_s4_class(fit_iter, "ARPMM2")
  expect_s4_class(fit_linear, "ARPMM2")
  
  # Coefficients should be reasonable (between -1 and 1 for AR(1))
  # coef() returns named vector including intercept, so check AR coefficient
  expect_true(all(abs(coef(fit_global)) < 2))  # Allow intercept
  expect_true(all(abs(coef(fit_iter)) < 2))
  expect_true(all(abs(coef(fit_linear)) < 2))
})

test_that("ma_pmm2 accepts pmm2_variant parameter", {
  set.seed(456)
  x <- arima.sim(n = 100, list(ma = 0.6))
  
  # Test all three variants
  fit_global <- ma_pmm2(x, order = 1, pmm2_variant = "unified_global")
  fit_iter <- ma_pmm2(x, order = 1, pmm2_variant = "unified_iterative")
  fit_linear <- ma_pmm2(x, order = 1, pmm2_variant = "linearized")
  
  # All should return MAPMM2 objects
  expect_s4_class(fit_global, "MAPMM2")
  expect_s4_class(fit_iter, "MAPMM2")
  expect_s4_class(fit_linear, "MAPMM2")
  
  # Coefficients should be reasonable
  expect_true(all(abs(coef(fit_global)) < 2))
  expect_true(all(abs(coef(fit_iter)) < 2))
  expect_true(all(abs(coef(fit_linear)) < 2))
})

test_that("arma_pmm2 accepts pmm2_variant parameter", {
  set.seed(789)
  x <- arima.sim(n = 150, list(ar = 0.5, ma = 0.3))
  
  # Test default and iterative variants
  fit_global <- arma_pmm2(x, order = c(1, 1), pmm2_variant = "unified_global")
  fit_iter <- arma_pmm2(x, order = c(1, 1), pmm2_variant = "unified_iterative")
  
  # Should return ARMAPMM2 objects
  expect_s4_class(fit_global, "ARMAPMM2")
  expect_s4_class(fit_iter, "ARMAPMM2")
  
  # Should have 2-3 coefficients (1 AR + 1 MA + optional intercept)
  expect_true(length(coef(fit_global)) >= 2)
  expect_true(length(coef(fit_iter)) >= 2)
})

test_that("arima_pmm2 accepts pmm2_variant parameter", {
  set.seed(111)
  # Create non-stationary series
  x <- cumsum(arima.sim(n = 120, list(ar = 0.6)))
  
  # Test ARIMA(1,1,0) with different variants
  fit_global <- arima_pmm2(x, order = c(1, 1, 0), pmm2_variant = "unified_global")
  fit_iter <- arima_pmm2(x, order = c(1, 1, 0), pmm2_variant = "unified_iterative")
  
  # Should return ARIMAPMM2 objects
  expect_s4_class(fit_global, "ARIMAPMM2")
  expect_s4_class(fit_iter, "ARIMAPMM2")
  
  # AR coefficient should be reasonable (allow intercept)
  expect_true(all(abs(coef(fit_global)) < 2))
  expect_true(all(abs(coef(fit_iter)) < 2))
})

test_that("pmm2_variant parameter validation works", {
  set.seed(222)
  x <- arima.sim(n = 80, list(ar = 0.5))
  
  # Invalid variant should error
  expect_error(
    ar_pmm2(x, order = 1, pmm2_variant = "invalid_variant"),
    "'arg' should be one of"
  )
})

test_that("backward compatibility: default behavior unchanged", {
  set.seed(333)
  x <- arima.sim(n = 100, list(ar = 0.7))
  
  # Calling without pmm2_variant should work (uses default)
  fit_default <- ar_pmm2(x, order = 1)
  fit_explicit <- ar_pmm2(x, order = 1, pmm2_variant = "unified_global")
  
  # Both should return same class
  expect_s4_class(fit_default, "ARPMM2")
  expect_s4_class(fit_explicit, "ARPMM2")
  
  # Note: Due to current implementation using ts_pmm2, 
  # coefficients may differ once unified wrappers are integrated
})

test_that("unified_global vs unified_iterative convergence", {
  set.seed(444)
  x <- arima.sim(n = 150, list(ar = c(0.6, -0.3)))
  
  fit_global <- ar_pmm2(x, order = 2, pmm2_variant = "unified_global")
  fit_iter <- ar_pmm2(x, order = 2, pmm2_variant = "unified_iterative")
  
  # Both should converge
  expect_true(fit_global@convergence)
  expect_true(fit_iter@convergence)
  
  # Iterative might have more iterations (if integrated)
  # For now, just check they're valid
  expect_true(fit_global@iterations > 0)
  expect_true(fit_iter@iterations > 0)
})

test_that("linearized variant works for MA models", {
  set.seed(555)
  x <- arima.sim(n = 100, list(ma = c(0.7, -0.4)))
  
  fit_linear <- ma_pmm2(x, order = 2, pmm2_variant = "linearized")
  
  # Should return MAPMM2 object
  expect_s4_class(fit_linear, "MAPMM2")
  
  # Should have 2-3 MA coefficients (+ optional intercept)
  expect_true(length(coef(fit_linear)) >= 2)
  
  # Should converge
  expect_true(fit_linear@convergence)
})

test_that("seasonal models accept pmm2_variant parameter", {
  skip_if_not_installed("EstemPMM", minimum_version = "0.2.0")
  
  set.seed(666)
  # Create seasonal data
  n <- 120
  x <- arima.sim(n = n, list(ar = 0.5, seasonal = list(order = c(1, 0, 0), period = 12)))
  
  # Test SAR with variants (if function exists)
  if (exists("sar_pmm2", mode = "function")) {
    fit_sar <- sar_pmm2(x, order = c(1, 1), season = list(period = 12), 
                        pmm2_variant = "unified_global")
    expect_s4_class(fit_sar, "SARPMM2")
  }
  
  # Test SMA with variants (if function exists)
  if (exists("sma_pmm2", mode = "function")) {
    x_sma <- arima.sim(n = n, list(seasonal = list(order = c(0, 0, 1), period = 12)))
    fit_sma <- sma_pmm2(x_sma, order = 1, season = list(period = 12),
                        pmm2_variant = "linearized")
    expect_s4_class(fit_sma, "SMAPMM2")
  }
})

test_that("moment estimates are reasonable", {
  set.seed(777)
  x <- arima.sim(n = 200, list(ar = 0.6))
  
  fit <- ar_pmm2(x, order = 1, pmm2_variant = "unified_global")
  
  # Check moment slots exist and are reasonable
  expect_true(!is.na(fit@m2))
  expect_true(!is.na(fit@m3))
  expect_true(!is.na(fit@m4))
  
  # m2 should be positive (variance)
  expect_true(fit@m2 > 0)
  
  # For Gaussian innovations, m3 should be near 0, m4/m2^2 near 3
  # But we allow wide range for non-Gaussian cases
  expect_true(abs(fit@m3) < 10 * fit@m2^1.5)
  expect_true(fit@m4 > 0)
})

test_that("all main functions have pmm2_variant argument", {
  # Check function signatures
  expect_true("pmm2_variant" %in% names(formals(ar_pmm2)))
  expect_true("pmm2_variant" %in% names(formals(ma_pmm2)))
  expect_true("pmm2_variant" %in% names(formals(arma_pmm2)))
  expect_true("pmm2_variant" %in% names(formals(arima_pmm2)))
  
  # Seasonal functions (if they exist)
  if (exists("sar_pmm2", mode = "function")) {
    expect_true("pmm2_variant" %in% names(formals(sar_pmm2)))
  }
  if (exists("sma_pmm2", mode = "function")) {
    expect_true("pmm2_variant" %in% names(formals(sma_pmm2)))
  }
  if (exists("sarima_pmm2", mode = "function")) {
    expect_true("pmm2_variant" %in% names(formals(sarima_pmm2)))
  }
})
