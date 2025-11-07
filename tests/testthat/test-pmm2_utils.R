set.seed(42)

test_that("pmm_skewness calculates expected shape", {
  x <- rnorm(100)
  skew <- pmm_skewness(x)

  expect_true(is.numeric(skew))
  expect_length(skew, 1)
  expect_false(is.na(skew))
  expect_true(abs(skew) < 0.5)
})

test_that("pmm_skewness distinguishes asymmetric data", {
  x_gamma <- rgamma(100, shape = 2)
  x_normal <- rnorm(100)

  expect_gt(pmm_skewness(x_gamma), pmm_skewness(x_normal))
})

test_that("pmm_kurtosis calculates expected value", {
  x <- rnorm(100)
  kurt <- pmm_kurtosis(x)

  expect_true(is.numeric(kurt))
  expect_length(kurt, 1)
  expect_false(is.na(kurt))
  expect_true(abs(kurt) < 1)
  kurt_with_baseline <- pmm_kurtosis(x, excess = FALSE)
  expect_true(abs(kurt_with_baseline - 3) < 1)
})

test_that("compute_moments returns valid moments", {
  x <- rnorm(100)
  moments <- compute_moments(x)

  expect_true(is.list(moments))
  expect_true(all(c("m2", "m3", "m4") %in% names(moments)))
  expect_true(moments$m2 > 0)
})

test_that("compute_moments handles constant series", {
  x_const <- rep(5, 100)
  moments <- compute_moments(x_const)

  expect_true(moments$m2 >= 0)
})
