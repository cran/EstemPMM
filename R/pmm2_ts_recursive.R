# pmm2_ts_recursive.R - Recursive (non-linearised) PMM2 estimation for models
# with a moving-average component.
#
# Background
# ----------
# The historical PMM2 time-series path builds the design matrix ONCE from the
# CSS residuals (see arma_build_design / ma_build_design) and keeps it fixed
# while the PMM2 Newton iteration runs.  For the AR columns this is exact: the
# regressors are observed data and do not depend on the parameter.  For the MA
# columns it is a ONE-STEP LINEARISATION of the innovation recursion, and it
# costs a factor of exactly (1 - theta^2) in asymptotic variance relative to
# CSS.  Under Gaussian innovations the linearised PMM2 estimator of an MA(1)
# parameter therefore has relative efficiency 1 - theta^2 < 1 against CSS, and
# the penalty does NOT vanish as the sample grows.
#
# The functions in this file recompute eps_t(beta) and the score regressors
# x_t(beta) = -d eps_t / d beta by the exact recursion at EVERY candidate
# parameter value.  When mu3 = 0 the PMM2 estimating equation then collapses
# identically onto the CSS first-order condition, so the (1 - theta^2) penalty
# disappears exactly rather than asymptotically.
#
# Model convention (matching stats::arima):
#   w_t = mu + sum_j phi_j w_{t-j} + eps_t + sum_k theta_k eps_{t-k}
#   eps_t = w_t - mu - sum_j phi_j w_{t-j} - sum_k theta_k eps_{t-k}
# with w_{s} = eps_{s} = 0 for s <= 0 (the CSS convention).


#' Exact ARMA innovation recursion with first and second derivatives
#'
#' Computes the innovation sequence \eqn{\varepsilon_t(\beta)} of an
#' ARMA(p, q) model together with the score regressors
#' \eqn{x_{t,r} = -\partial \varepsilon_t / \partial \beta_r} and, optionally,
#' the second-derivative array
#' \eqn{G_{t,rs} = \partial x_{t,r} / \partial \beta_s
#'      = -\partial^2 \varepsilon_t / \partial \beta_r \partial \beta_s}.
#'
#' @details
#' All three quantities obey the same linear recursion driven by different
#' inputs, so each is obtained with one call to \code{stats::filter}:
#' \deqn{\varepsilon_t = d_t - \sum_m \theta_m \varepsilon_{t-m},\quad
#'       d_t = w_t - \mu - \sum_j \phi_j w_{t-j}}
#' \deqn{x_{t,r} = b_{t,r} - \sum_m \theta_m x_{t-m,r}}
#' with \eqn{b_{t,r} = 1} for the intercept, \eqn{w_{t-j}} for \eqn{\phi_j} and
#' \eqn{\varepsilon_{t-k}} for \eqn{\theta_k}, and
#' \deqn{G_{t,rs} = -[r=\theta_k]\,x_{t-k,s} - [s=\theta_l]\,x_{t-l,r}
#'                 - \sum_m \theta_m G_{t-m,rs}.}
#' Note that for \eqn{q = 0} the recursion is the identity and
#' \eqn{x_{t,\phi_j} = w_{t-j}}, i.e. the recursive construction reduces
#' exactly to the fixed design matrix used by the linearised path.  The two
#' solvers therefore agree bit-for-bit on pure AR models.
#'
#' @param w Numeric vector: the (differenced, if applicable) series.
#' @param phi Numeric vector of AR coefficients (length p, possibly empty).
#' @param theta Numeric vector of MA coefficients (length q, possibly empty).
#' @param mu Numeric intercept (only used when \code{include_intercept = TRUE}).
#' @param include_intercept Logical: treat \code{mu} as a free parameter and
#'   return its column in \code{X}.
#' @param second_order Logical: also return the second-derivative array
#'   \code{G}.
#'
#' @return List with \code{eps} (length n), \code{X} (n x npar), \code{G}
#'   (n x npar x npar array or \code{NULL}) and \code{par_names}.
#'   Parameter order is \code{c(intercept, phi_1..phi_p, theta_1..theta_q)},
#'   the intercept being present only when requested.
#' @keywords internal
arma_recursion <- function(w, phi = numeric(0), theta = numeric(0),
                           mu = 0, include_intercept = FALSE,
                           second_order = FALSE) {
  w <- as.numeric(w)
  n <- length(w)
  p <- length(phi)
  q <- length(theta)
  k_int <- if (isTRUE(include_intercept)) 1L else 0L
  npar <- k_int + p + q

  if (npar == 0L) {
    stop("arma_recursion(): model has no free parameters")
  }

  # y_t = v_t - sum_m theta_m y_{t-m}, with zero initial conditions
  rec <- function(v) {
    if (q == 0L) {
      return(as.numeric(v))
    }
    as.numeric(stats::filter(v, filter = -theta, method = "recursive"))
  }
  lagv <- function(v, k) {
    if (k <= 0L) v else c(numeric(k), v[seq_len(max(n - k, 0L))])
  }

  drive <- w
  if (k_int == 1L) drive <- drive - mu
  if (p > 0L) {
    for (j in seq_len(p)) drive <- drive - phi[j] * lagv(w, j)
  }
  eps <- rec(drive)

  X <- matrix(0, nrow = n, ncol = npar)
  if (k_int == 1L) X[, 1L] <- rec(rep(1, n))
  if (p > 0L) {
    for (j in seq_len(p)) X[, k_int + j] <- rec(lagv(w, j))
  }
  if (q > 0L) {
    for (k in seq_len(q)) X[, k_int + p + k] <- rec(lagv(eps, k))
  }

  G <- NULL
  if (isTRUE(second_order)) {
    G <- array(0, dim = c(n, npar, npar))
    if (q > 0L) {
      # index of theta_k within the parameter vector
      is_theta <- function(r) r > k_int + p
      theta_lag <- function(r) r - k_int - p
      for (r in seq_len(npar)) {
        for (s in seq_len(r)) {
          rhs <- numeric(n)
          if (is_theta(r)) rhs <- rhs - lagv(X[, s], theta_lag(r))
          if (is_theta(s)) rhs <- rhs - lagv(X[, r], theta_lag(s))
          if (any(rhs != 0)) {
            val <- rec(rhs)
            G[, r, s] <- val
            if (s != r) G[, s, r] <- val
          }
        }
      }
    }
  }

  par_names <- c(
    if (k_int == 1L) "intercept" else NULL,
    if (p > 0L) paste0("ar", seq_len(p)) else NULL,
    if (q > 0L) paste0("ma", seq_len(q)) else NULL
  )

  list(eps = eps, X = X, G = G, par_names = par_names,
       p = p, q = q, k_int = k_int)
}


#' PMM2 estimating function and Jacobian on the exact recursion
#'
#' @param par Parameter vector \code{c(intercept, phi, theta)}.
#' @param w Series (differenced/centred as appropriate).
#' @param p,q AR and MA orders.
#' @param include_intercept Logical.
#' @param m2,m3,m4 Central moments of the innovations (held fixed).
#' @param n_burn Number of leading observations excluded from the sums.
#' @param jacobian Logical: also return the Jacobian.
#' @param exact_jacobian Logical: include the second-derivative term.  When
#'   \code{FALSE} the Gauss-Newton approximation is used.
#'
#' @return List with \code{g} (npar), \code{J} (npar x npar or \code{NULL}) and
#'   \code{eps}.  The sign convention matches \code{\link{pmm2_algorithm}}:
#'   the Newton step is \code{par - solve(J, g)}.
#' @keywords internal
pmm2_recursive_score <- function(par, w, p, q, include_intercept,
                                 m2, m3, m4, n_burn = p,
                                 jacobian = TRUE, exact_jacobian = TRUE) {
  k_int <- if (isTRUE(include_intercept)) 1L else 0L
  mu <- if (k_int == 1L) par[1L] else 0
  phi <- if (p > 0L) par[k_int + seq_len(p)] else numeric(0)
  theta <- if (q > 0L) par[k_int + p + seq_len(q)] else numeric(0)

  rc <- arma_recursion(w, phi, theta, mu = mu,
                       include_intercept = include_intercept,
                       second_order = jacobian && exact_jacobian && q > 0L)

  n <- length(w)
  keep <- if (n_burn > 0L) seq.int(n_burn + 1L, n) else seq_len(n)
  eps <- rc$eps[keep]
  X <- rc$X[keep, , drop = FALSE]

  A <- m4 - m2^2
  # Same polynomial as pmm2_algorithm(), rewritten in terms of eps = y - X b:
  #   Z1 = m3 eps^2 - (m4 - m2^2) eps - m2 m3
  Z1 <- m3 * eps^2 - A * eps - m2 * m3
  g <- as.numeric(crossprod(X, Z1))

  J <- NULL
  if (isTRUE(jacobian)) {
    wjac <- A - 2 * m3 * eps
    J <- crossprod(X, X * wjac)
    if (exact_jacobian && !is.null(rc$G)) {
      Gk <- rc$G[keep, , , drop = FALSE]
      npar <- ncol(X)
      for (r in seq_len(npar)) {
        for (s in seq_len(r)) {
          add <- sum(Z1 * Gk[, r, s])
          J[r, s] <- J[r, s] + add
          if (s != r) J[s, r] <- J[s, r] + add
        }
      }
    }
  }

  list(g = g, J = J, eps = rc$eps, X_full = rc$X)
}


#' Check ARMA stationarity / invertibility
#'
#' @param phi AR coefficients.
#' @param theta MA coefficients.
#' @param tol Root modulus must exceed \code{1 + tol}.
#' @return \code{TRUE} if all roots lie outside the unit circle.
#' @keywords internal
arma_admissible <- function(phi = numeric(0), theta = numeric(0), tol = 1e-6) {
  ok_poly <- function(coefs) {
    if (length(coefs) == 0L) return(TRUE)
    if (any(!is.finite(coefs))) return(FALSE)
    # polynomial 1 + coefs[1] z + ... ; roots must be outside the unit circle
    rr <- tryCatch(polyroot(c(1, coefs)), error = function(e) NULL)
    if (is.null(rr)) return(FALSE)
    all(Mod(rr) > 1 + tol)
  }
  ok_poly(-phi) && ok_poly(theta)
}


#' Solve the recursive PMM2 estimating equations
#'
#' Safeguarded Newton iteration on the PMM2 estimating equations built from the
#' exact ARMA innovation recursion.  Globalisation combines
#' \itemize{
#'   \item the analytic Jacobian (second recursion), with a Gauss-Newton
#'         fallback when the exact Jacobian is singular or produces no descent;
#'   \item backtracking line search on \eqn{\|g\|_2};
#'   \item projection onto the stationary / invertible region: any step that
#'         leaves it is halved.
#' }
#' If no admissible decreasing step exists the initial (CSS-based) value is
#' returned with \code{convergence = FALSE}, so the solver never returns
#' \code{NA}.
#'
#' @param w Series.
#' @param par_init Starting value \code{c(intercept, phi, theta)}.
#' @param p,q AR and MA orders.
#' @param include_intercept Logical.
#' @param m2,m3,m4 Central moments of the innovations (held fixed).
#' @param max_iter Maximum Newton iterations.
#' @param tol Convergence tolerance on the parameter step.
#' @param n_burn Leading observations excluded from the estimating sums.
#'   Defaults to \code{p}, matching \code{stats::arima(method = "CSS")}.
#' @param max_halvings Maximum backtracking halvings per iteration.
#' @param verbose Print progress.
#'
#' @return List with \code{par}, \code{eps}, \code{convergence},
#'   \code{iterations}, \code{gnorm} and \code{fallback}.
#' @keywords internal
pmm2_recursive_solve <- function(w, par_init, p, q, include_intercept,
                                 m2, m3, m4,
                                 max_iter = 50, tol = 1e-8,
                                 n_burn = p, max_halvings = 25L,
                                 verbose = FALSE) {
  k_int <- if (isTRUE(include_intercept)) 1L else 0L
  split_par <- function(par) {
    list(
      phi = if (p > 0L) par[k_int + seq_len(p)] else numeric(0),
      theta = if (q > 0L) par[k_int + p + seq_len(q)] else numeric(0)
    )
  }
  admissible <- function(par) {
    if (any(!is.finite(par))) return(FALSE)
    sp <- split_par(par)
    arma_admissible(sp$phi, sp$theta)
  }
  score <- function(par, jac) {
    pmm2_recursive_score(par, w, p, q, include_intercept,
                         m2, m3, m4, n_burn = n_burn,
                         jacobian = jac, exact_jacobian = TRUE)
  }

  par <- as.numeric(par_init)
  if (!admissible(par)) {
    # Nudge an inadmissible start back inside the region.
    sp <- split_par(par)
    shrink <- 1
    for (i in seq_len(40L)) {
      shrink <- shrink * 0.9
      cand <- par
      if (p > 0L) cand[k_int + seq_len(p)] <- sp$phi * shrink
      if (q > 0L) cand[k_int + p + seq_len(q)] <- sp$theta * shrink
      if (admissible(cand)) {
        par <- cand
        break
      }
    }
  }

  cur <- score(par, jac = TRUE)
  if (any(!is.finite(cur$g))) {
    return(list(par = as.numeric(par_init), eps = NULL, convergence = FALSE,
                iterations = 0L, gnorm = NA_real_, fallback = TRUE))
  }
  gnorm <- sqrt(sum(cur$g^2))

  converged <- FALSE
  iterations <- 0L
  fallback <- FALSE

  for (iter in seq_len(max_iter)) {
    iterations <- iter

    step <- tryCatch(solve(cur$J, cur$g), error = function(e) NULL)
    used_gn <- FALSE
    if (is.null(step) || any(!is.finite(step))) {
      # Gauss-Newton fallback: drop the second-derivative term.
      gn <- pmm2_recursive_score(par, w, p, q, include_intercept,
                                 m2, m3, m4, n_burn = n_burn,
                                 jacobian = TRUE, exact_jacobian = FALSE)
      step <- tryCatch(solve(gn$J, gn$g), error = function(e) NULL)
      used_gn <- TRUE
    }
    if (is.null(step) || any(!is.finite(step))) {
      fallback <- TRUE
      break
    }

    alpha <- 1
    accepted <- FALSE
    for (h in seq_len(max_halvings)) {
      cand <- par - alpha * step
      if (admissible(cand)) {
        trial <- score(cand, jac = FALSE)
        if (all(is.finite(trial$g))) {
          gnew <- sqrt(sum(trial$g^2))
          if (gnew < gnorm * (1 - 1e-4 * alpha)) {
            accepted <- TRUE
            break
          }
        }
      }
      alpha <- alpha / 2
    }

    if (!accepted) {
      # Try the Gauss-Newton direction once before giving up.
      if (!used_gn) {
        gn <- pmm2_recursive_score(par, w, p, q, include_intercept,
                                   m2, m3, m4, n_burn = n_burn,
                                   jacobian = TRUE, exact_jacobian = FALSE)
        step2 <- tryCatch(solve(gn$J, gn$g), error = function(e) NULL)
        if (!is.null(step2) && all(is.finite(step2))) {
          alpha <- 1
          for (h in seq_len(max_halvings)) {
            cand <- par - alpha * step2
            if (admissible(cand)) {
              trial <- score(cand, jac = FALSE)
              if (all(is.finite(trial$g))) {
                gnew <- sqrt(sum(trial$g^2))
                if (gnew < gnorm * (1 - 1e-4 * alpha)) {
                  step <- step2
                  accepted <- TRUE
                  break
                }
              }
            }
            alpha <- alpha / 2
          }
        }
      }
    }

    if (!accepted) {
      # Stationary point of ||g|| that is not a root, or numerically flat:
      # stop and keep the best parameter found so far.
      converged <- gnorm < 1e-6 * max(1, length(w))
      break
    }

    delta <- alpha * step
    par <- par - delta
    cur <- score(par, jac = TRUE)
    gnorm <- sqrt(sum(cur$g^2))

    if (verbose) {
      cat("  iter", iter, "| step =", formatC(sqrt(sum(delta^2)), digits = 6),
          "| ||g|| =", formatC(gnorm, digits = 6), "\n")
    }

    if (sqrt(sum(delta^2)) < tol) {
      converged <- TRUE
      break
    }
  }

  if (any(!is.finite(par))) {
    par <- as.numeric(par_init)
    converged <- FALSE
    fallback <- TRUE
  }

  sp <- split_par(par)
  mu <- if (k_int == 1L) par[1L] else 0
  eps <- arma_recursion(w, sp$phi, sp$theta, mu = mu,
                        include_intercept = include_intercept)$eps

  list(par = par, eps = eps, convergence = converged,
       iterations = iterations, gnorm = gnorm, fallback = fallback)
}


#' Score regressors of a fitted recursive PMM2 time-series model
#'
#' Rebuilds the exact score regressors
#' \eqn{x_t = -\partial \varepsilon_t / \partial \beta} at the fitted parameter
#' value. These are the \eqn{\mathcal F_{t-1}}-measurable regressors that enter
#' the closed-form sandwich covariance; unlike the frozen CSS design matrix used
#' by the linearised path, they are the true derivatives of the criterion.
#'
#' @param object A \code{TS2fit} object of model type \code{"ma"},
#'   \code{"arma"} or \code{"arima"}.
#'
#' @return List with \code{X} (effective rows x number of slope parameters) and
#'   \code{n_eff}.
#' @keywords internal
pmm2_recursive_design <- function(object) {
  mt <- object@model_type
  if (!mt %in% c("ar", "ma", "arma", "arima")) {
    stop("pmm2_recursive_design() supports non-seasonal AR/MA/ARMA/ARIMA fits only")
  }
  ord <- object@order
  getord <- function(nm) if (is.null(ord[[nm]])) 0L else as.integer(ord[[nm]])
  p <- getord("ar")
  q <- getord("ma")
  d <- getord("d")

  x <- as.numeric(object@original_series)
  w <- if (d > 0L) diff(x, differences = d) else x
  w <- w - as.numeric(object@intercept)

  cf <- as.numeric(object@coefficients)
  phi <- if (p > 0L) cf[seq_len(p)] else numeric(0)
  theta <- if (q > 0L) cf[p + seq_len(q)] else numeric(0)

  rc <- arma_recursion(w, phi, theta, include_intercept = FALSE)
  n <- length(w)
  keep <- if (p > 0L) seq.int(p + 1L, n) else seq_len(n)
  X <- rc$X[keep, , drop = FALSE]
  colnames(X) <- c(
    if (p > 0L) paste0("ar", seq_len(p)) else NULL,
    if (q > 0L) paste0("ma", seq_len(q)) else NULL
  )
  list(X = X, n_eff = length(keep))
}

#' Recursive PMM2 fit for a pure MA(q) model
#'
#' Drop-in replacement for the linearised \code{ma_pmm2_fit()} that recomputes
#' the innovation recursion at every candidate parameter value.
#'
#' @param x Original series.
#' @param q MA order.
#' @param css_fit Initial CSS fit (from \code{ma_css_fit}).
#' @param include.mean Logical: estimate an intercept correction.
#' @param max_iter,tol Solver controls.
#' @param verbose Print progress.
#'
#' @return List with \code{coefficients}, \code{intercept}, \code{innovations},
#'   \code{convergence} and \code{iterations}.
#' @keywords internal
ma_pmm2_fit_recursive <- function(x, q, css_fit, include.mean = TRUE,
                                  max_iter = 50, tol = 1e-8, verbose = FALSE) {
  x <- as.numeric(x)
  moments <- compute_moments(css_fit$residuals)
  w <- x - css_fit$intercept

  k_int <- if (isTRUE(include.mean)) 1L else 0L
  par_init <- c(if (k_int == 1L) 0 else NULL, css_fit$coefficients)

  res <- pmm2_recursive_solve(
    w = w, par_init = par_init, p = 0L, q = q,
    include_intercept = isTRUE(include.mean),
    m2 = moments$m2, m3 = moments$m3, m4 = moments$m4,
    max_iter = max_iter, tol = tol, n_burn = 0L, verbose = verbose
  )

  mu_hat <- if (k_int == 1L) res$par[1L] else 0
  theta <- res$par[k_int + seq_len(q)]
  intercept <- css_fit$intercept + mu_hat
  innovations <- ma_compute_innovations(x - intercept, theta, q)

  list(
    coefficients = theta,
    intercept = intercept,
    innovations = innovations,
    convergence = res$convergence,
    iterations = res$iterations
  )
}
