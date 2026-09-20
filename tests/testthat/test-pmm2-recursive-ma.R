# Tests for the recursive (non-linearised) treatment of the MA recursion.
# See R/pmm2_ts_recursive.R and NEWS.md (0.5.0).

# ---- helpers ---------------------------------------------------------------

sim_ma1 <- function(n, theta, gen = rnorm, m = 0, s = 1, burn = 200L) {
  e <- (gen(n + burn) - m) / s
  w <- e + theta * c(0, utils::head(e, -1))
  w[(burn + 1L):(n + burn)]
}

css_objective <- function(par, w, p, q, n_burn) {
  e <- EstemPMM:::arma_recursion(
    w,
    phi   = if (p > 0) par[seq_len(p)] else numeric(0),
    theta = if (q > 0) par[p + seq_len(q)] else numeric(0)
  )$eps
  sum(e[(n_burn + 1L):length(e)]^2)
}

# ---- (0) the recursion itself ----------------------------------------------

test_that("arma_recursion reproduces the analytic derivatives", {
  skip_if_not_installed("numDeriv")
  set.seed(1)
  w <- rnorm(60)
  par0 <- c(0.3, 0.4, -0.2, 0.5, -0.3) # mu, phi1, phi2, theta1, theta2
  rc <- EstemPMM:::arma_recursion(w, phi = par0[2:3], theta = par0[4:5],
                                  mu = par0[1], include_intercept = TRUE,
                                  second_order = TRUE)
  f_eps <- function(par) {
    EstemPMM:::arma_recursion(w, phi = par[2:3], theta = par[4:5], mu = par[1],
                              include_intercept = TRUE)$eps
  }
  # x_t = -d eps_t / d beta
  expect_lt(max(abs(rc$X + numDeriv::jacobian(f_eps, par0))), 1e-7)

  # G_{t,r,s} = d x_{t,r} / d beta_s
  worst <- 0
  for (r in seq_len(5)) {
    f_X <- function(par) {
      EstemPMM:::arma_recursion(w, phi = par[2:3], theta = par[4:5], mu = par[1],
                                include_intercept = TRUE)$X[, r]
    }
    worst <- max(worst, max(abs(numDeriv::jacobian(f_X, par0) - rc$G[, r, ])))
  }
  expect_lt(worst, 1e-6)
})

test_that("with q = 0 the recursion degenerates to the fixed lag design", {
  set.seed(2)
  w <- rnorm(50)
  rc <- EstemPMM:::arma_recursion(w, phi = c(0.5, -0.2), theta = numeric(0))
  lagm <- cbind(c(0, w[1:49]), c(0, 0, w[1:48]))
  expect_equal(rc$X, lagm, ignore_attr = TRUE)
})

# ---- (a) exact identity with CSS under symmetric innovations ---------------
#
# When mu3 = 0 the PMM2 estimating equation is
#   sum_t x_t(theta) eps_t(theta) = 0,
# which is EXACTLY the CSS first-order condition.  This is an identity, not an
# approximation, so the tolerance here is numerical, not statistical.

test_that("recursive PMM2 with m3 = 0 is identical to CSS (MA(1))", {
  set.seed(7)
  worst <- 0
  for (rep in seq_len(10)) {
    w <- sim_ma1(600, theta = -0.5)
    fit <- stats::arima(w, order = c(0, 0, 1), include.mean = FALSE,
                        method = "CSS")
    mm <- EstemPMM:::compute_moments(as.numeric(fit$residuals))
    sol <- EstemPMM:::pmm2_recursive_solve(
      w, par_init = 0, p = 0L, q = 1L, include_intercept = FALSE,
      m2 = mm$m2, m3 = 0, m4 = mm$m4, max_iter = 60L, tol = 1e-12, n_burn = 0L
    )
    expect_true(sol$convergence)
    # against the same CSS criterion, minimised to machine precision
    obj <- function(par) css_objective(par, w = w, p = 0L, q = 1L, n_burn = 0L)
    opt <- stats::optim(as.numeric(fit$coef[1]), obj, method = "BFGS",
                        control = list(reltol = 1e-14))
    worst <- max(worst, abs(sol$par - opt$par))
  }
  expect_lt(worst, 1e-6)
})

test_that("recursive PMM2 with m3 = 0 is identical to CSS (ARMA(1,1))", {
  set.seed(11)
  worst <- 0
  for (rep in seq_len(6)) {
    x <- as.numeric(stats::arima.sim(n = 800, list(ar = 0.6, ma = -0.4)))
    fit <- stats::arima(x, order = c(1, 0, 1), include.mean = FALSE,
                        method = "CSS")
    start <- as.numeric(fit$coef[c("ar1", "ma1")])
    mm <- EstemPMM:::compute_moments(as.numeric(fit$residuals))
    sol <- EstemPMM:::pmm2_recursive_solve(
      x, par_init = start, p = 1L, q = 1L, include_intercept = FALSE,
      m2 = mm$m2, m3 = 0, m4 = mm$m4, max_iter = 60L, tol = 1e-12, n_burn = 1L
    )
    expect_true(sol$convergence)
    obj <- function(par) css_objective(par, w = x, p = 1L, q = 1L, n_burn = 1L)
    opt <- stats::optim(start, obj, method = "BFGS",
                        control = list(reltol = 1e-14))
    worst <- max(worst, max(abs(sol$par - opt$par)))
  }
  expect_lt(worst, 1e-5)
})

# ---- (b) regression test: the linearised branch keeps its 1 - theta^2 -------

test_that("linearised branch is exactly OLS on the frozen design when m3 = 0", {
  # Structural form of the penalty: with m3 = 0 the fixed-design PMM2 solves
  #   sum_t xhat_t (w_t - xhat_t' b) = 0,
  # i.e. plain OLS of w on the FROZEN CSS residual lags.  It is this one-step
  # linearisation of the recursion that costs the factor 1 - theta^2.
  set.seed(21)
  w <- sim_ma1(500, theta = -0.5)
  fit <- stats::arima(w, order = c(0, 0, 1), include.mean = FALSE,
                      method = "CSS")
  e_css <- EstemPMM:::arma_recursion(w, theta = as.numeric(fit$coef[1]))$eps
  d <- EstemPMM:::arma_build_design(w, e_css, p = 0, q = 1, intercept = 0,
                                    include_intercept = FALSE)
  algo <- EstemPMM:::pmm2_algorithm(
    b_init = as.numeric(fit$coef[1]), X = d$X, y = d$y,
    m2 = stats::var(e_css), m3 = 0, m4 = mean(e_css^4),
    max_iter = 50L, tol = 1e-12, regularize = FALSE, reg_lambda = 0
  )
  ols <- as.numeric(stats::coef(stats::lm(d$y ~ d$X - 1)))
  expect_equal(algo$b, ols, tolerance = 1e-8)
})

test_that("linearised branch pays the 1 - theta^2 penalty, recursive does not", {
  skip_on_cran()
  set.seed(4242)
  theta0 <- -0.5
  reps <- 400L
  est <- matrix(NA_real_, reps, 3,
                dimnames = list(NULL, c("css", "lin", "rec")))
  for (r in seq_len(reps)) {
    w <- sim_ma1(500, theta = theta0)
    fit <- tryCatch(stats::arima(w, order = c(0, 0, 1), include.mean = FALSE,
                                 method = "CSS"), error = function(e) NULL)
    if (is.null(fit)) next
    th <- as.numeric(fit$coef[1])
    est[r, "css"] <- th
    e_css <- EstemPMM:::arma_recursion(w, theta = th)$eps
    mm <- EstemPMM:::compute_moments(e_css)
    d <- EstemPMM:::arma_build_design(w, e_css, p = 0, q = 1, intercept = 0,
                                      include_intercept = FALSE)
    a <- tryCatch(EstemPMM:::pmm2_algorithm(
      b_init = th, X = d$X, y = d$y, m2 = mm$m2, m3 = mm$m3, m4 = mm$m4,
      max_iter = 50L, tol = 1e-8, regularize = FALSE, reg_lambda = 0),
      error = function(e) NULL)
    if (!is.null(a)) est[r, "lin"] <- a$b[1]
    est[r, "rec"] <- EstemPMM:::pmm2_recursive_solve(
      w, par_init = th, p = 0L, q = 1L, include_intercept = FALSE,
      m2 = mm$m2, m3 = mm$m3, m4 = mm$m4, max_iter = 50L, tol = 1e-10,
      n_burn = 0L)$par
  }
  v <- apply(est, 2, function(z) stats::var(z[is.finite(z)]))
  re_lin <- v[["css"]] / v[["lin"]]
  re_rec <- v[["css"]] / v[["rec"]]
  # Frozen expectation for the linearised branch under Gaussian innovations.
  expect_gt(re_lin / (1 - theta0^2), 0.80)
  expect_lt(re_lin / (1 - theta0^2), 1.25)
  # The recursive branch removes the penalty.
  expect_gt(re_rec, 0.90)
  expect_gt(re_rec, re_lin)
})

# ---- (c) solver convergence ------------------------------------------------

test_that("recursive solver converges on every replication of the MA grid", {
  skip_on_cran()
  set.seed(12345)
  laws <- list(
    Gaussian      = list(gen = rnorm, m = 0, s = 1),
    Gamma         = list(gen = function(n) rgamma(n, 2, 1), m = 2, s = sqrt(2)),
    `Chi-squared` = list(gen = function(n) rchisq(n, 3), m = 3, s = sqrt(6))
  )
  n_try <- 0L; n_finite <- 0L; n_conv <- 0L
  for (L in laws) {
    for (theta0 in c(-0.8, -0.5, -0.2, 0.4)) {
      for (r in seq_len(25L)) {
        w <- sim_ma1(300, theta = theta0, gen = L$gen, m = L$m, s = L$s)
        fit <- tryCatch(stats::arima(w, order = c(0, 0, 1),
                                     include.mean = FALSE, method = "CSS"),
                        error = function(e) NULL)
        if (is.null(fit)) next
        th <- as.numeric(fit$coef[1])
        mm <- EstemPMM:::compute_moments(
          EstemPMM:::arma_recursion(w, theta = th)$eps)
        sol <- EstemPMM:::pmm2_recursive_solve(
          w, par_init = th, p = 0L, q = 1L, include_intercept = FALSE,
          m2 = mm$m2, m3 = mm$m3, m4 = mm$m4, max_iter = 50L, tol = 1e-10,
          n_burn = 0L)
        n_try <- n_try + 1L
        n_finite <- n_finite + is.finite(sol$par)
        n_conv <- n_conv + isTRUE(sol$convergence)
        # never leaves the invertible region
        expect_true(EstemPMM:::arma_admissible(theta = sol$par))
      }
    }
  }
  expect_gt(n_try, 250L)
  expect_equal(n_finite, n_try)
  expect_equal(n_conv, n_try)
})

test_that("recursive solver never returns NA from an inadmissible start", {
  set.seed(99)
  w <- sim_ma1(300, theta = 0.7)
  for (start in c(-0.99, -0.5, 0, 0.5, 0.99, 1.5, -2)) {
    mm <- EstemPMM:::compute_moments(EstemPMM:::arma_recursion(w, theta = 0.7)$eps)
    sol <- EstemPMM:::pmm2_recursive_solve(
      w, par_init = start, p = 0L, q = 1L, include_intercept = FALSE,
      m2 = mm$m2, m3 = mm$m3, m4 = mm$m4, max_iter = 50L, tol = 1e-10,
      n_burn = 0L)
    expect_true(is.finite(sol$par))
    expect_true(EstemPMM:::arma_admissible(theta = sol$par))
  }
})

# ---- integration through the user-facing API -------------------------------

test_that("ma_solver is plumbed through and recorded on the fit", {
  set.seed(31)
  w <- sim_ma1(300, theta = 0.6, gen = function(n) rgamma(n, 2, 1),
               m = 2, s = sqrt(2))
  f_lin <- ma_pmm2(w, order = 1)
  f_rec <- ma_pmm2(w, order = 1, ma_solver = "recursive")
  expect_identical(f_lin@ma_solver, "linearized")
  expect_identical(f_rec@ma_solver, "recursive")
  expect_false(isTRUE(all.equal(coef(f_lin), coef(f_rec))))
})

test_that("the default remains the linearised branch", {
  set.seed(32)
  w <- sim_ma1(200, theta = 0.5)
  expect_identical(ma_pmm2(w, order = 1)@ma_solver, "linearized")
  expect_identical(arma_pmm2(w, order = c(1, 1))@ma_solver, "linearized")
  expect_identical(arima_pmm2(cumsum(w), order = c(0, 1, 1))@ma_solver,
                   "linearized")
})

test_that("pure AR fits are bit-identical under both settings", {
  set.seed(33)
  x <- as.numeric(stats::arima.sim(n = 300, list(ar = c(0.5, -0.2))))
  a1 <- ts_pmm2(x, order = 2, model_type = "ar")
  a2 <- ts_pmm2(x, order = 2, model_type = "ar", ma_solver = "recursive")
  expect_identical(coef(a1), coef(a2))
  b1 <- arima_pmm2(cumsum(x), order = c(2, 1, 0))
  b2 <- arima_pmm2(cumsum(x), order = c(2, 1, 0), ma_solver = "recursive")
  expect_identical(coef(b1), coef(b2))
})

test_that("recursive branch handles q > 1 and mixed ARIMA(p, d, q)", {
  set.seed(34)
  x2 <- as.numeric(stats::arima.sim(n = 600, list(ma = c(0.5, -0.3))))
  f2 <- ma_pmm2(x2, order = 2, ma_solver = "recursive")
  expect_length(f2@coefficients, 2L)
  expect_true(all(is.finite(f2@coefficients)))
  expect_true(EstemPMM:::arma_admissible(theta = f2@coefficients))

  x3 <- as.numeric(stats::arima.sim(n = 800,
                                    list(ar = c(0.5, -0.3), ma = c(0.4, 0.2))))
  f3 <- arma_pmm2(x3, order = c(2, 2), ma_solver = "recursive")
  expect_length(f3@coefficients, 4L)
  expect_true(all(is.finite(f3@coefficients)))

  f4 <- arima_pmm2(cumsum(x3), order = c(1, 1, 1), ma_solver = "recursive")
  expect_length(f4@coefficients, 2L)
  expect_true(all(is.finite(f4@coefficients)))
})

# ---- closed-form sandwich covariance ---------------------------------------

test_that("vcov()/confint() work for recursive MA fits and refuse otherwise", {
  set.seed(35)
  w <- sim_ma1(600, theta = 0.6, gen = function(n) rgamma(n, 2, 1),
               m = 2, s = sqrt(2))
  f_rec <- ma_pmm2(w, order = 1, ma_solver = "recursive", include.mean = FALSE)
  V <- vcov(f_rec)
  expect_equal(dim(V), c(1L, 1L))
  expect_true(is.finite(V[1, 1]) && V[1, 1] > 0)
  expect_identical(colnames(V), "ma1")
  ci <- confint(f_rec)
  expect_true(ci[1, 1] < coef(f_rec)[1] && coef(f_rec)[1] < ci[1, 2])

  f_lin <- ma_pmm2(w, order = 1, include.mean = FALSE)
  expect_error(vcov(f_lin), "ma_solver")
  expect_error(confint(f_lin), "ma_solver")
})

test_that("closed-form sandwich agrees with the numerical Jacobian sandwich", {
  set.seed(36)
  w <- sim_ma1(1000, theta = -0.5, gen = function(n) rgamma(n, 2, 1),
               m = 2, s = sqrt(2))
  f <- ma_pmm2(w, order = 1, ma_solver = "recursive", include.mean = FALSE)
  V_closed <- vcov(f)[1, 1]

  # Numerical sandwich A^{-1} B A^{-T} built from the empirical psi_t.
  th <- as.numeric(coef(f))
  rc <- EstemPMM:::arma_recursion(w, theta = th, second_order = TRUE)
  a <- f@m4 - f@m2^2
  eps <- rc$eps
  x <- rc$X[, 1]
  B_t <- a * eps - f@m3 * (eps^2 - f@m2)
  psi <- x * B_t
  n <- length(w)
  A_hat <- mean(rc$G[, 1, 1] * B_t - x^2 * (a - 2 * f@m3 * eps))
  B_hat <- mean(psi^2)
  V_num <- B_hat / (A_hat^2) / n
  expect_equal(V_closed, V_num, tolerance = 0.05)
})
