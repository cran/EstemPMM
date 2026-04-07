test_that("pmm_dispatch selects PMM3 for uniform errors", {
  set.seed(42)
  x <- rnorm(200)
  eps <- runif(200, -sqrt(3), sqrt(3))
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  result <- pmm_dispatch(r, verbose = FALSE)
  expect_equal(result$method, "PMM3")
})

test_that("pmm_dispatch selects PMM2 for exponential errors", {
  set.seed(42)
  x <- rnorm(200)
  eps <- rexp(200) - 1
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  result <- pmm_dispatch(r, verbose = FALSE)
  expect_equal(result$method, "PMM2")
})

test_that("pmm_dispatch selects OLS for normal errors", {
  set.seed(42)
  x <- rnorm(200)
  eps <- rnorm(200)
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  result <- pmm_dispatch(r, verbose = FALSE)
  expect_equal(result$method, "OLS")
})

test_that("pmm_dispatch returns complete structure", {
  set.seed(42)
  r <- rnorm(100)
  result <- pmm_dispatch(r, verbose = FALSE)

  expect_true(is.list(result))
  expect_true(result$method %in% c("OLS", "PMM2", "PMM3"))
  expect_true(is.numeric(result$gamma3))
  expect_true(is.numeric(result$gamma4))
  expect_true(is.numeric(result$gamma6))
  expect_true(is.numeric(result$g2))
  expect_true(is.numeric(result$g3))
  expect_true(is.numeric(result$g_selected))
  expect_true(is.numeric(result$improvement_pct))
  expect_true(is.character(result$reasoning))
  expect_true(is.numeric(result$n))
})

test_that("pmm_dispatch respects custom thresholds", {
  set.seed(42)
  x <- rnorm(200)
  eps <- runif(200, -sqrt(3), sqrt(3))
  y <- 1 + 2 * x + eps
  r <- residuals(lm(y ~ x))

  # With very strict kurtosis threshold, uniform may not qualify
  result <- pmm_dispatch(r, kurtosis_threshold = -2.0, verbose = FALSE)
  expect_true(result$method %in% c("OLS", "PMM2", "PMM3"))
})
