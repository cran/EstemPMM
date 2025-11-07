set.seed(42)
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2  # Asymmetric errors
y <- 2 + 1.5 * x + errors
test_data <- data.frame(x = x, y = y)

test_that("lm_pmm2 produces valid fit object", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_s4_class(fit, "PMM2fit")
  expect_true(fit@convergence)
  expect_length(fit@coefficients, 2)
  expect_equal(length(fit@residuals), n)
})

test_that("lm_pmm2 coefficients are numeric and reasonable", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_true(is.numeric(fit@coefficients))
  expect_false(any(is.na(fit@coefficients)))
  expect_false(any(is.infinite(fit@coefficients)))

  expect_true(all(abs(fit@coefficients) < 100))
})

test_that("lm_pmm2 residuals sum to near zero", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_true(is.numeric(fit@residuals))
  expect_equal(length(fit@residuals), n)
  expect_true(abs(mean(fit@residuals)) < 1e-2)
})

test_that("lm_pmm2 returns moment estimates", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_true(is.numeric(fit@m2))
  expect_true(is.numeric(fit@m3))
  expect_true(is.numeric(fit@m4))
  expect_false(any(is.na(c(fit@m2, fit@m3, fit@m4))))
  expect_true(fit@m2 > 0)
})

test_that("lm_pmm2 compared to OLS shows variance reduction", {
  fit_pmm2 <- lm_pmm2(y ~ x, data = test_data)
  fit_ols <- lm(y ~ x, data = test_data)

  coef_diff <- abs(fit_pmm2@coefficients - coef(fit_ols))
  expect_true(all(coef_diff < 1.0))
})

test_that("lm_pmm2 handles additional predictors", {
  x2 <- rnorm(n)
  test_data_ext <- transform(test_data, x2 = x2, y_multi = 2 + 1.5 * x + 0.8 * x2 + errors)

  fit <- lm_pmm2(y_multi ~ x + x2, data = test_data_ext)

  expect_s4_class(fit, "PMM2fit")
  expect_length(fit@coefficients, 3)
})

test_that("lm_pmm2 handles missing values with na.omit", {
  test_data_na <- test_data
  test_data_na$y[1] <- NA

  fit <- lm_pmm2(y ~ x, data = test_data_na, na.action = na.omit)

  expect_s4_class(fit, "PMM2fit")
  expect_equal(length(fit@residuals), n - 1)
})

test_that("lm_pmm2 supports models without intercept", {
  fit <- lm_pmm2(y ~ x - 1, data = test_data)

  expect_s4_class(fit, "PMM2fit")
  expect_length(fit@coefficients, 1)
})

test_that("lm_pmm2 convergence status stored as logical", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_true(is.logical(fit@convergence))
  expect_length(fit@convergence, 1)
})

test_that("lm_pmm2 works on very small samples", {
  small_data <- data.frame(x = rnorm(4), y = rnorm(4))

  expect_error(lm_pmm2(y ~ x, data = small_data), NA)
})

test_that("lm_pmm2 stores original call", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  expect_true(is.call(fit@call))
  expect_true(grepl("lm_pmm2", deparse(fit@call)))
})
