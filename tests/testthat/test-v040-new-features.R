# test-v040-new-features.R
#
# Tests for features introduced in EstemPMM 0.4.0 in response to the
# JSS reviewer report:
#   * Unified API: pmm_lm / pmm_ar / pmm_ma / pmm_arma / pmm_arima /
#                  pmm_sarima
#   * S3 class PMMdispatch: print / summary / format
#   * show() methods for PMM2fit / PMM3fit / TS2fit / TS3fit
#   * ci_method argument and bias column in pmm2_inference /
#     ts_pmm2_inference

## ---- shared fixtures --------------------------------------------------------
set.seed(2026L)
n   <- 150
x   <- rnorm(n)
eps <- rgamma(n, shape = 2, scale = 1) - 2   # skewed: gamma3 > 0.3
y   <- 1.5 + 2.0 * x + eps
dat <- data.frame(y = y, x = x)

eps_sym <- runif(n, -2, 2)                    # symmetric platykurtic
y_sym   <- 0.5 + 1.5 * x + eps_sym
dat_sym <- data.frame(y = y_sym, x = x)

x_ts <- as.numeric(arima.sim(list(ar = 0.6), n = 200,
                              rand.gen = function(k) rgamma(k, 2, 1) - 2))

# =============================================================================
# 1. Unified regression API: pmm_lm()
# =============================================================================

test_that("pmm_lm(method='pmm2') returns PMM2fit", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  expect_s4_class(fit, "PMM2fit")
  expect_true(is(fit, "PMMfit"))
  expect_true(fit@convergence)
})

test_that("pmm_lm(method='pmm3') returns PMM3fit", {
  fit <- suppressWarnings(pmm_lm(y ~ x, data = dat_sym, method = "pmm3"))
  expect_s4_class(fit, "PMM3fit")
  expect_true(is(fit, "PMMfit"))
})

test_that("pmm_lm(method='auto') dispatches to PMM2 for skewed data", {
  fit <- pmm_lm(y ~ x, data = dat, method = "auto")
  # Gamma residuals have |gamma3| > 0.3 so dispatcher should pick PMM2
  expect_s4_class(fit, "PMM2fit")
})

test_that("pmm_lm(method='auto') dispatches to PMM3 for platykurtic data", {
  fit <- suppressWarnings(pmm_lm(y ~ x, data = dat_sym, method = "auto"))
  # Uniform residuals: gamma4 ~ -1.2 < -0.7 and gamma3 ~ 0; dispatch -> PMM3
  expect_s4_class(fit, "PMM3fit")
})

test_that("pmm_lm preserves the outer call so show() and predict() work", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  # The call slot must reference pmm_lm, not lm_pmm2, so that
  # predict() can eval(call$formula) correctly.
  call_fn <- as.character(fit@call[[1L]])
  expect_equal(call_fn, "pmm_lm")
})

test_that("pmm_lm result supports all generic extractors", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  expect_true(is.numeric(coef(fit)))
  expect_true(is.numeric(residuals(fit)))
  expect_true(is.numeric(fitted(fit)))
  expect_equal(nobs(fit), n)
  expect_s3_class(logLik(fit), "logLik")
  expect_true(is.finite(AIC(fit)))
  expect_true(is.finite(BIC(fit)))
})

# =============================================================================
# 2. Unified time-series API: pmm_ar / pmm_ma / pmm_arma / pmm_arima
# =============================================================================

test_that("pmm_ar(method='pmm2') returns ARPMM2 inheriting PMMtsfit", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm2")
  expect_s4_class(fit, "ARPMM2")
  expect_true(is(fit, "PMMtsfit"))
  expect_true(is(fit, "PMMfit"))
  expect_true(fit@convergence)
})

test_that("pmm_ar(method='pmm3') returns ARPMM3 inheriting PMMtsfit", {
  fit <- ar_pmm3(x_ts, order = 1)         # call via legacy name too
  fit2 <- pmm_ar(x_ts, order = 1, method = "pmm3")
  expect_s4_class(fit2, "ARPMM3")
  expect_true(is(fit2, "PMMtsfit"))
})

test_that("pmm_arima(method='pmm2') returns ARIMAPMM2", {
  x_arima <- cumsum(x_ts)
  fit <- pmm_arima(x_arima, order = c(1, 1, 0), method = "pmm2")
  expect_s4_class(fit, "ARIMAPMM2")
  expect_true(is(fit, "TS2fit"))
})

test_that("pmm_arima(method='pmm3') returns ARIMAPMM3", {
  x_arima <- cumsum(x_ts)
  fit <- suppressWarnings(pmm_arima(x_arima, order = c(1, 1, 0), method = "pmm3"))
  expect_s4_class(fit, "ARIMAPMM3")
  expect_true(is(fit, "TS3fit"))
})

test_that("pmm_sarima with unsupported method='pmm3' errors informatively", {
  expect_error(pmm_sarima(x_ts, order = c(1, 0, 0),
                          seasonal = list(order = c(1, 0, 0), period = 4),
                          method = "pmm3"),
               regexp = "not implemented")
})

test_that("pmm_ar preserves call for show() output", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm2")
  expect_equal(as.character(fit@call[[1L]]), "pmm_ar")
})

# =============================================================================
# 3. S3 class PMMdispatch: print / summary / format
# =============================================================================

test_that("pmm_dispatch returns PMMdispatch S3 object", {
  res <- pmm_dispatch(residuals(lm(y ~ x, data = dat)), verbose = FALSE)
  expect_s3_class(res, "PMMdispatch")
})

test_that("PMMdispatch $ accessor returns expected fields", {
  res <- pmm_dispatch(residuals(lm(y ~ x, data = dat)), verbose = FALSE)
  expect_true(res$method %in% c("OLS", "PMM2", "PMM3"))
  expect_true(is.numeric(res$gamma3))
  expect_true(is.numeric(res$gamma4))
  expect_true(is.numeric(res$g2))
  expect_true(is.numeric(res$g3))
  expect_true(is.numeric(res$improvement_pct))
  expect_true(is.character(res$reasoning))
  expect_true(is.integer(res$n) || is.numeric(res$n))
})

test_that("print.PMMdispatch produces human-readable output", {
  res <- pmm_dispatch(residuals(lm(y ~ x, data = dat)), verbose = FALSE)
  out <- capture.output(print(res))
  expect_true(any(grepl("dispatch", out, ignore.case = TRUE)))
  expect_true(any(grepl("gamma3", out, ignore.case = TRUE)))
  expect_true(any(grepl("Selected", out, ignore.case = TRUE)))
})

test_that("summary.PMMdispatch returns a single-row data frame", {
  res <- pmm_dispatch(residuals(lm(y ~ x, data = dat)), verbose = FALSE)
  s   <- summary(res)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 1L)
  expect_true("method" %in% names(s))
  expect_true("g_selected" %in% names(s))
  expect_true("improvement_pct" %in% names(s))
})

test_that("format.PMMdispatch returns a single character string", {
  res <- pmm_dispatch(residuals(lm(y ~ x, data = dat)), verbose = FALSE)
  f   <- format(res)
  expect_type(f, "character")
  expect_length(f, 1L)
  expect_true(grepl("PMMdispatch", f))
})

test_that("pmm_dispatch verbose=TRUE calls print.PMMdispatch (no error)", {
  expect_output(pmm_dispatch(residuals(lm(y ~ x, data = dat)), verbose = TRUE),
                regexp = "dispatch|gamma|Selected",
                ignore.case = TRUE)
})

# =============================================================================
# 4. show() methods for fit classes
# =============================================================================

test_that("show(PMM2fit) prints lm-style header, not slot dump", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  out <- capture.output(show(fit))
  # Must contain "Call:" and "Coefficients:", must NOT start with "An object"
  expect_true(any(grepl("Call:", out, fixed = TRUE)))
  expect_true(any(grepl("Coefficients:", out, fixed = TRUE)))
  expect_false(any(grepl("^An object of class", out)))
  # Must not show raw slot names like @coefficients / Slot "coefficients"
  expect_false(any(grepl("Slot", out, fixed = TRUE)))
})

test_that("show(PMM2fit) footer contains convergence information", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  out <- paste(capture.output(show(fit)), collapse = " ")
  expect_true(grepl("iteration|converge", out, ignore.case = TRUE))
})

test_that("show(PMM3fit) works without slot dump", {
  fit <- suppressWarnings(pmm_lm(y ~ x, data = dat_sym, method = "pmm3"))
  out <- capture.output(show(fit))
  expect_true(any(grepl("Call:", out, fixed = TRUE)))
  expect_true(any(grepl("Coefficients:", out, fixed = TRUE)))
  expect_false(any(grepl("^An object of class", out)))
})

test_that("show(TS2fit) displays model label, not slot dump", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm2")
  out <- capture.output(show(fit))
  # Must mention the model type
  expect_true(any(grepl("AR", out, ignore.case = TRUE)))
  expect_false(any(grepl("^An object of class", out)))
})

test_that("show(TS3fit) displays model label", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm3")
  out <- capture.output(show(fit))
  expect_true(any(grepl("AR|Model", out, ignore.case = TRUE)))
  expect_false(any(grepl("^An object of class", out)))
})

test_that("typing fit at console triggers show(), not default slot dump", {
  # Verify that the class has a show method so the REPL uses it
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  has_show <- isVirtualClass("PMM2fit") ||
              !is.null(getMethod("show", "PMM2fit", optional = TRUE))
  expect_true(has_show)
})

# =============================================================================
# 5. ci_method argument and bias column in pmm2_inference
# =============================================================================

test_that("pmm2_inference returns bias column", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  out <- pmm2_inference(fit, y ~ x, data = dat, B = 50, seed = 42L)
  expect_true("bias" %in% names(out))
  expect_true(is.numeric(out$bias))
})

test_that("pmm2_inference ci_method='normal' CI contains the estimate", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  out <- pmm2_inference(fit, y ~ x, data = dat, B = 80, seed = 42L,
                        ci_method = "normal")
  # Normal CI is symmetric about the point estimate — must always contain it
  expect_true(all(out$conf.low  <= out$Estimate + 1e-12))
  expect_true(all(out$conf.high >= out$Estimate - 1e-12))
})

test_that("pmm2_inference ci_method='percentile' is the old default", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  out_p <- pmm2_inference(fit, y ~ x, data = dat, B = 80, seed = 42L,
                          ci_method = "percentile")
  out_n <- pmm2_inference(fit, y ~ x, data = dat, B = 80, seed = 42L,
                          ci_method = "normal")
  # The two methods should produce different interval endpoints
  expect_false(isTRUE(all.equal(out_p$conf.low, out_n$conf.low)))
})

test_that("pmm2_inference ci_method='basic' pivotal interval is reflected", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  out_b <- pmm2_inference(fit, y ~ x, data = dat, B = 80, seed = 42L,
                          ci_method = "basic")
  out_p <- pmm2_inference(fit, y ~ x, data = dat, B = 80, seed = 42L,
                          ci_method = "percentile")
  # basic:   [2*est - q97.5, 2*est - q2.5]
  # percentile: [q2.5, q97.5]
  # => basic_low + percentile_high ≈ 2 * est (within MC noise)
  expect_equal(out_b$conf.low + out_p$conf.high,
               2 * out_b$Estimate,
               tolerance = 0.01)
})

test_that("pmm2_inference ci_method defaults to 'normal'", {
  fit <- pmm_lm(y ~ x, data = dat, method = "pmm2")
  out_default <- pmm2_inference(fit, y ~ x, data = dat, B = 80, seed = 42L)
  out_normal  <- pmm2_inference(fit, y ~ x, data = dat, B = 80, seed = 42L,
                                ci_method = "normal")
  expect_equal(out_default$conf.low,  out_normal$conf.low)
  expect_equal(out_default$conf.high, out_normal$conf.high)
})

# =============================================================================
# 6. ci_method and bias in ts_pmm2_inference
# =============================================================================

test_that("ts_pmm2_inference returns bias column", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm2")
  out <- ts_pmm2_inference(fit, B = 50, method = "residual", seed = 42L)
  expect_true("bias" %in% names(out))
  expect_true(is.numeric(out$bias))
})

test_that("ts_pmm2_inference normal CI contains the estimate (residual boot)", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm2")
  out <- ts_pmm2_inference(fit, B = 80, method = "residual", seed = 42L,
                           ci_method = "normal")
  expect_true(all(out$conf.low  <= out$Estimate + 1e-12))
  expect_true(all(out$conf.high >= out$Estimate - 1e-12))
})

test_that("ts_pmm2_inference normal CI contains the estimate (block boot)", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm2")
  out <- ts_pmm2_inference(fit, B = 80, method = "block",
                           block_length = 7L, seed = 42L,
                           ci_method = "normal")
  # Normal CI is Wald-based (est ± z * SE); must always contain est
  expect_true(all(out$conf.low  <= out$Estimate + 1e-12))
  expect_true(all(out$conf.high >= out$Estimate - 1e-12))
})

test_that("ts_pmm2_inference ci_method defaults to 'normal'", {
  fit <- pmm_ar(x_ts, order = 1, method = "pmm2")
  out_d <- ts_pmm2_inference(fit, B = 60, method = "residual", seed = 99L)
  out_n <- ts_pmm2_inference(fit, B = 60, method = "residual", seed = 99L,
                             ci_method = "normal")
  expect_equal(out_d$conf.low,  out_n$conf.low)
  expect_equal(out_d$conf.high, out_n$conf.high)
})

test_that("ts_pmm2_inference block-boot percentile CI may miss estimate (known bias)", {
  # This is a known finite-sample property for block bootstrap on
  # AR coefficients near the stationarity boundary; the test documents
  # that the behaviour is acknowledged, not accidentally broken.
  set.seed(123L)
  x_persistent <- arima.sim(list(ar = 0.85), n = 300,
                            rand.gen = function(k) rgamma(k, 2, 1) - 2)
  fit <- pmm_ar(x_persistent, order = 1, method = "pmm2")
  out <- ts_pmm2_inference(fit, B = 100, method = "block",
                           block_length = 10L, seed = 42L,
                           ci_method = "percentile")
  # Row-names may vary; use column-vector indexing robustly.
  bias_val  <- out$bias[[1L]]
  cl_val    <- out$conf.low[[1L]]
  ch_val    <- out$conf.high[[1L]]
  # The bias should be visibly negative (boot distribution shrinks toward 0).
  expect_true(bias_val < -0.01)
  # The percentile CI may or may not contain the estimate — we just check
  # that bias and CI columns are present and finite.
  expect_true(is.finite(cl_val))
  expect_true(is.finite(ch_val))
})
