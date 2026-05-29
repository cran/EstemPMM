# test-s4-generics.R
# Tests for the four new S4 generics: vcov, confint, logLik, nobs

set.seed(42)
n   <- 150
x   <- rnorm(n)
eps <- rgamma(n, shape = 2, scale = 1) - 2   # skewed, gamma3 ~ 1.41
y   <- 1.5 + 2.8 * x + eps
dat <- data.frame(y = y, x = x)

fit_lm   <- lm_pmm2(y ~ x, data = dat)
fit_lm3  <- lm_pmm3(y ~ x, data = dat)

# ---- logLik ------------------------------------------------------------------

test_that("logLik.PMM2fit returns logLik object", {
  ll <- logLik(fit_lm)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_equal(attr(ll, "nobs"), n)
  expect_equal(attr(ll, "df"), length(coef(fit_lm)))
})

test_that("logLik.PMM2fit is consistent with AIC (AIC = -2*ll + 2*p)", {
  ll  <- logLik(fit_lm)
  aic <- AIC(fit_lm)
  p   <- length(coef(fit_lm))
  expect_equal(aic, -2 * as.numeric(ll) + 2 * p, tolerance = 1e-8)
})

# ---- nobs --------------------------------------------------------------------

test_that("nobs.PMM2fit returns integer equal to sample size", {
  expect_equal(nobs(fit_lm), n)
})

test_that("nobs.PMM2fit equals length of residuals", {
  expect_equal(nobs(fit_lm), length(residuals(fit_lm)))
})

# ---- vcov --------------------------------------------------------------------

test_that("vcov.PMM2fit returns a matrix", {
  V <- vcov(fit_lm)
  expect_true(is.matrix(V))
  expect_equal(nrow(V), 2L)
  expect_equal(ncol(V), 2L)
})

test_that("vcov.PMM2fit is symmetric and positive-definite", {
  V <- vcov(fit_lm)
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(diag(V) > 0))
  expect_true(det(V) > 0)
})

test_that("vcov.PMM2fit is smaller than OLS vcov under skewed errors", {
  V_pmm2 <- vcov(fit_lm)
  V_ols  <- vcov(lm(y ~ x, data = dat))
  # For gamma errors g2 < 1, so PMM2 variance < OLS variance
  expect_true(V_pmm2[2, 2] <= V_ols[2, 2] * 1.1)  # allow 10% tolerance
})

test_that("vcov.PMM2fit has correct row/col names", {
  V   <- vcov(fit_lm)
  nms <- names(coef(fit_lm))
  expect_equal(rownames(V), nms)
  expect_equal(colnames(V), nms)
  expect_true(all(nms %in% c("(Intercept)", "x")))
})

test_that("vcov.PMM2fit errors when model_matrix is missing", {
  fit_bare <- new("PMM2fit",
                  coefficients = c(1, 2),
                  residuals    = rnorm(10),
                  m2 = 1, m3 = 0.5, m4 = 3,
                  convergence  = TRUE,
                  iterations   = 5L,
                  call         = quote(lm_pmm2(y ~ x)))
  expect_error(vcov(fit_bare), "model_matrix not found")
})

# ---- confint -----------------------------------------------------------------

test_that("confint.PMM2fit returns a matrix with correct dimensions", {
  ci <- confint(fit_lm)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2L)
  expect_equal(ncol(ci), 2L)
})

test_that("confint.PMM2fit lower < upper", {
  ci <- confint(fit_lm)
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("confint.PMM2fit 95% CI contains true values (approx)", {
  # True beta1 = 2.8; with n=150 this CI should cover it
  ci <- confint(fit_lm)
  expect_true(ci["x", 1] < 2.8 && 2.8 < ci["x", 2])
})

test_that("confint.PMM2fit level argument changes interval width", {
  ci99 <- confint(fit_lm, level = 0.99)
  ci90 <- confint(fit_lm, level = 0.90)
  width99 <- ci99["x", 2] - ci99["x", 1]
  width90 <- ci90["x", 2] - ci90["x", 1]
  expect_true(width99 > width90)
})

test_that("confint.PMM2fit parm argument subsets output", {
  ci <- confint(fit_lm, parm = "x")
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "x")
})

# ---- TS2fit (AR model) -------------------------------------------------------

set.seed(42)
ar_series <- arima.sim(list(ar = 0.7), n = 300,
                       rand.gen = function(n) rgamma(n, 2, 1) - 2)
fit_ar  <- ar_pmm2(ar_series, order = 1)
fit_ar2 <- ar_pmm2(ar_series, order = 2)

test_that("logLik.TS2fit returns logLik object", {
  ll <- logLik(fit_ar)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
  expect_equal(attr(ll, "nobs"), nobs(fit_ar))
})

test_that("nobs.TS2fit returns effective sample size", {
  expect_equal(nobs(fit_ar), length(residuals(fit_ar)))
  expect_true(nobs(fit_ar) <= 300L)
})

test_that("vcov.TS2fit for AR(1) returns 1x1 matrix", {
  V <- vcov(fit_ar)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(1L, 1L))
  expect_true(V[1, 1] > 0)
  expect_equal(rownames(V), "ar1")
})

test_that("vcov.TS2fit for AR(2) returns 2x2 matrix", {
  V <- vcov(fit_ar2)
  expect_equal(dim(V), c(2L, 2L))
  expect_true(all(diag(V) > 0))
})

test_that("confint.TS2fit for AR(1) contains true value", {
  ci <- confint(fit_ar)
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "ar1")
  expect_true(ci["ar1", 1] < 0.7 && 0.7 < ci["ar1", 2])
})

test_that("vcov.TS2fit errors for non-AR models", {
  fit_arima <- arima_pmm2(ar_series, order = c(1, 1, 0))
  expect_error(vcov(fit_arima), "AR models")
})

test_that("confint.TS2fit errors for non-AR models with helpful message", {
  fit_ma <- ma_pmm2(ar_series, order = 1)
  expect_error(confint(fit_ma), "ts_pmm2_inference")
})

# ---- PMM3fit vcov/confint ----------------------------------------------------

# Symmetric platykurtic errors: uniform on [-a, a] has gamma4 = -1.2
set.seed(2026)
n_sym   <- 300
x_sym   <- rnorm(n_sym)
eps_sym <- runif(n_sym, -2, 2)            # symmetric, gamma4 < 0
y_sym   <- 0.5 + 1.5 * x_sym + eps_sym
dat_sym <- data.frame(y = y_sym, x = x_sym)

fit_p3 <- suppressWarnings(lm_pmm3(y ~ x, data = dat_sym))

test_that("vcov.PMM3fit returns a matrix with correct dimensions", {
  V <- vcov(fit_p3)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(2L, 2L))
})

test_that("vcov.PMM3fit is symmetric and positive-definite", {
  V <- vcov(fit_p3)
  expect_equal(V, t(V), tolerance = 1e-12)
  expect_true(all(diag(V) > 0))
  expect_true(det(V) > 0)
})

test_that("vcov.PMM3fit is smaller than OLS vcov under platykurtic errors", {
  V_p3  <- vcov(fit_p3)
  V_ols <- vcov(lm(y ~ x, data = dat_sym))
  # For uniform errors gamma4 = -1.2, so g3 < 1 and PMM3 variance < OLS variance
  expect_true(V_p3[2, 2] < V_ols[2, 2])
})

test_that("vcov.PMM3fit equals g3 * vcov(OLS) under same residuals", {
  V_p3  <- vcov(fit_p3)
  V_ols <- vcov(lm(y ~ x, data = dat_sym))
  # The PMM3 efficiency factor (computed from initial OLS residuals)
  # must reproduce the OLS-PMM3 ratio on the diagonal.
  g3 <- fit_p3@g_coefficient
  expect_equal(V_p3[2, 2] / V_ols[2, 2], g3, tolerance = 0.05)
})

test_that("vcov.PMM3fit has correct row/col names", {
  V   <- vcov(fit_p3)
  nms <- names(coef(fit_p3))
  expect_equal(rownames(V), nms)
  expect_equal(colnames(V), nms)
})

test_that("vcov.PMM3fit errors when model_matrix is missing", {
  fit_bare <- new("PMM3fit",
                  coefficients  = c(1, 2),
                  residuals     = rnorm(10),
                  m2 = 1, m4 = 3, m6 = 15,
                  gamma4 = 0, gamma6 = 0,
                  g_coefficient = 1, kappa = NA_real_,
                  convergence  = TRUE,
                  iterations   = 5L,
                  call         = quote(lm_pmm3(y ~ x)))
  expect_error(vcov(fit_bare), "model_matrix not found")
})

test_that("confint.PMM3fit returns a matrix with correct dimensions", {
  ci <- confint(fit_p3)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2L)
  expect_equal(ncol(ci), 2L)
  expect_true(all(ci[, 1] < ci[, 2]))
})

test_that("confint.PMM3fit 95% CI contains true value", {
  ci <- confint(fit_p3)
  expect_true(ci["x", 1] < 1.5 && 1.5 < ci["x", 2])
})

test_that("confint.PMM3fit level argument changes interval width", {
  ci99 <- confint(fit_p3, level = 0.99)
  ci90 <- confint(fit_p3, level = 0.90)
  expect_true((ci99["x", 2] - ci99["x", 1]) > (ci90["x", 2] - ci90["x", 1]))
})

# ---- TS3fit vcov/confint (AR only) -------------------------------------------

set.seed(2027)
ar_sym <- arima.sim(list(ar = 0.6), n = 400,
                    rand.gen = function(k) runif(k, -2, 2))
fit_ar3 <- ar_pmm3(ar_sym, order = 1)

test_that("vcov.TS3fit for AR(1) returns 1x1 positive matrix", {
  V <- vcov(fit_ar3)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(1L, 1L))
  expect_true(V[1, 1] > 0)
  expect_equal(rownames(V), "ar1")
})

test_that("confint.TS3fit for AR(1) contains true value", {
  ci <- confint(fit_ar3)
  expect_equal(nrow(ci), 1L)
  expect_equal(rownames(ci), "ar1")
  expect_true(ci["ar1", 1] < 0.6 && 0.6 < ci["ar1", 2])
})

test_that("vcov.TS3fit errors for non-AR PMM3 models", {
  fit_ma3 <- suppressWarnings(ma_pmm3(ar_sym, order = 1))
  expect_error(vcov(fit_ma3), "AR models")
})
