test_that("compute_moments_pmm3 returns correct structure for uniform data", {
  set.seed(42)
  r <- runif(500, -sqrt(3), sqrt(3))
  mom <- compute_moments_pmm3(r)

  expect_true(is.list(mom))
  expect_true(mom$m2 > 0)
  expect_true(mom$m4 > 0)
  expect_true(mom$m6 > 0)
  expect_true(mom$gamma4 < 0)        # Uniform is platykurtic
  expect_true(abs(mom$gamma3) < 0.3)  # Uniform is symmetric
  expect_true(mom$g3 < 1)             # Should show improvement
  expect_false(is.na(mom$kappa))
  expect_true(mom$improvement_pct > 0)
})

test_that("compute_moments_pmm3 handles normal data", {
  set.seed(42)
  r <- rnorm(500)
  mom <- compute_moments_pmm3(r)

  expect_true(abs(mom$gamma4) < 1)  # gamma4 near 0 for normal
  expect_true(abs(mom$gamma3) < 0.3)
})

test_that("pmm3_variance_factor returns correct structure", {
  vf <- pmm3_variance_factor(1.0, 1.8, 9.0)
  expect_true(is.list(vf))
  expect_true(!is.null(vf$gamma4))
  expect_true(!is.null(vf$gamma6))
  expect_true(!is.null(vf$g3))
})

test_that("pmm3_variance_factor handles m2 = 0", {
  vf <- pmm3_variance_factor(0, 0, 0)
  expect_true(is.na(vf$gamma4))
  expect_true(is.na(vf$gamma6))
  expect_true(is.na(vf$g3))
})

test_that("pmm3_variance_factor handles near-Gaussian", {
  # For standard normal: m2=1, m4=3, m6=15 => gamma4=0, gamma6=0
  vf <- pmm3_variance_factor(1, 3, 15)
  expect_equal(vf$gamma4, 0)
  expect_equal(vf$gamma6, 0)
  expect_equal(vf$g3, 1)  # No improvement for Gaussian
})

test_that("pmm_gamma6 returns near-zero for normal data", {
  set.seed(42)
  x <- rnorm(1000)
  g6 <- pmm_gamma6(x)
  expect_true(abs(g6) < 5)  # Approximately 0 for large samples
})

test_that("pmm_gamma6 returns positive for uniform data", {
  set.seed(42)
  x <- runif(1000, -1, 1)
  g6 <- pmm_gamma6(x)
  # Theoretical gamma6 for uniform ~ 6.86
  expect_true(g6 > 0)
})

test_that("pmm_gamma6 warns for too few values", {
  expect_warning(pmm_gamma6(1:5), "At least 6")
})

test_that("test_symmetry identifies symmetric data", {
  set.seed(42)
  x <- rnorm(200)
  result <- test_symmetry(x)
  expect_true(result$is_symmetric)
  expect_true(abs(result$gamma3) <= 0.3)
})

test_that("test_symmetry identifies asymmetric data", {
  set.seed(42)
  x <- rexp(200)
  result <- test_symmetry(x)
  expect_false(result$is_symmetric)
  expect_true(abs(result$gamma3) > 0.3)
})

test_that("test_symmetry respects custom threshold", {
  set.seed(42)
  x <- rnorm(200)
  result <- test_symmetry(x, threshold = 0.01)
  # With threshold = 0.01 even slight asymmetry from sampling triggers FALSE
  expect_true(is.logical(result$is_symmetric))
})
