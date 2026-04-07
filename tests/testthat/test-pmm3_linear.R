set.seed(42)
n <- 200
x <- rnorm(n)
errors_unif <- runif(n, -sqrt(3), sqrt(3))  # Symmetric platykurtic (gamma4 ~ -1.2)
y <- 2 + 1.5 * x + errors_unif
test_data <- data.frame(x = x, y = y)

test_that("lm_pmm3 produces valid fit object", {
  fit <- lm_pmm3(y ~ x, data = test_data)

  expect_s4_class(fit, "PMM3fit")
  expect_true(fit@convergence)
  expect_length(fit@coefficients, 2)
  expect_equal(length(fit@residuals), n)
})

test_that("lm_pmm3 coefficients are numeric and reasonable", {
  fit <- lm_pmm3(y ~ x, data = test_data)

  expect_true(is.numeric(fit@coefficients))
  expect_false(any(is.na(fit@coefficients)))
  expect_false(any(is.infinite(fit@coefficients)))
  expect_true(all(abs(fit@coefficients) < 100))
})

test_that("lm_pmm3 residuals sum to near zero", {
  fit <- lm_pmm3(y ~ x, data = test_data)

  expect_true(is.numeric(fit@residuals))
  expect_equal(length(fit@residuals), n)
  expect_true(abs(mean(fit@residuals)) < 0.1)
})

test_that("lm_pmm3 returns correct moment slots", {
  fit <- lm_pmm3(y ~ x, data = test_data)

  expect_true(fit@m2 > 0)
  expect_true(fit@m4 > 0)
  expect_true(fit@m6 > 0)
  expect_true(fit@gamma4 < 0)  # Uniform errors are platykurtic
  expect_true(fit@g_coefficient < 1)  # Should show improvement
})

test_that("lm_pmm3 with Gaussian errors converges and approximates OLS", {
  set.seed(123)
  y_gauss <- 2 + 1.5 * x + rnorm(n)
  gauss_data <- data.frame(x = x, y = y_gauss)

  fit <- lm_pmm3(y ~ x, data = gauss_data)
  fit_ols <- lm(y ~ x, data = gauss_data)

  expect_true(fit@convergence)
  expect_true(all(abs(coef(fit) - coef(fit_ols)) < 0.1))
})

test_that("lm_pmm3 warns for asymmetric errors", {
  set.seed(99)
  y_asym <- 2 + 1.5 * x + (rexp(n) - 1)
  asym_data <- data.frame(x = x, y = y_asym)

  expect_warning(lm_pmm3(y ~ x, data = asym_data), "asymmetric|PMM2")
})

test_that("lm_pmm3 handles multiple predictors", {
  x2 <- rnorm(n)
  y_multi <- 2 + 1.5 * x + 0.8 * x2 + errors_unif
  multi_data <- data.frame(x = x, x2 = x2, y = y_multi)

  fit <- lm_pmm3(y ~ x + x2, data = multi_data)

  expect_s4_class(fit, "PMM3fit")
  expect_length(fit@coefficients, 3)
})

test_that("lm_pmm3 handles missing values with na.omit", {
  test_data_na <- test_data
  test_data_na$y[1] <- NA

  fit <- lm_pmm3(y ~ x, data = test_data_na, na.action = na.omit)

  expect_s4_class(fit, "PMM3fit")
  expect_equal(length(fit@residuals), n - 1)
})

test_that("lm_pmm3 supports models without intercept", {
  fit <- lm_pmm3(y ~ x - 1, data = test_data)

  expect_s4_class(fit, "PMM3fit")
  expect_length(fit@coefficients, 1)
})

test_that("lm_pmm3 stores original call", {
  fit <- lm_pmm3(y ~ x, data = test_data)

  expect_true(is.call(fit@call))
  expect_true(grepl("lm_pmm3", deparse(fit@call)))
})

test_that("lm_pmm3 adaptive mode converges", {
  fit <- lm_pmm3(y ~ x, data = test_data, adaptive = TRUE)

  expect_s4_class(fit, "PMM3fit")
  expect_true(fit@convergence)
})

test_that("coef method works for PMM3fit", {
  fit <- lm_pmm3(y ~ x, data = test_data)
  expect_equal(coef(fit), fit@coefficients)
})

test_that("residuals method works for PMM3fit", {
  fit <- lm_pmm3(y ~ x, data = test_data)
  expect_equal(residuals(fit), fit@residuals)
})

test_that("fitted method works for PMM3fit", {
  fit <- lm_pmm3(y ~ x, data = test_data)
  fv <- fitted(fit)
  expect_equal(length(fv), n)
  expect_true(is.numeric(fv))
})

test_that("predict method works for PMM3fit with newdata", {
  fit <- lm_pmm3(y ~ x, data = test_data)
  nd <- data.frame(x = c(0, 1, 2))
  preds <- predict(fit, newdata = nd)
  expect_length(preds, 3)
  expect_true(is.numeric(preds))
})

test_that("summary method runs without error for PMM3fit", {
  fit <- lm_pmm3(y ~ x, data = test_data)
  expect_output(summary(fit), "PMM3 estimation results")
})
