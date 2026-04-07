# pmm3_solver.R - Internal Newton-Raphson solver for PMM3

#' PMM3 Newton-Raphson solver (vectorised)
#'
#' Solves the PMM3 score equations using Newton-Raphson iteration.
#' The score is \eqn{Z = X'(\varepsilon(\kappa - \varepsilon^2))} and the
#' Jacobian is \eqn{J = X' \mathrm{diag}(3\varepsilon^2 - \kappa) X}.
#'
#' @param B_ols     OLS starting values (length P)
#' @param X         Design matrix (n x P)
#' @param Y         Response vector (length n)
#' @param kappa     Moment ratio (scalar)
#' @param adaptive  Re-estimate kappa each iteration? (default FALSE)
#' @param tol       Convergence tolerance on ||delta||_2 (default 1e-6)
#' @param max_iter  Maximum NR iterations (default 100)
#' @param step_max  Maximum step size (default 5.0)
#'
#' @return List: B (estimates), converged, iter, kappa, near_gaussian
#' @keywords internal
.pmm3_nr_solver <- function(B_ols, X, Y, kappa,
                             adaptive = FALSE,
                             tol      = 1e-6,
                             max_iter = 100,
                             step_max = 5.0) {

  P <- length(B_ols)
  B <- B_ols

  # Near-Gaussian: kappa undefined -> return OLS
  if (is.na(kappa)) {
    return(list(B = B_ols, converged = TRUE, iter = 0L,
                kappa = NA_real_, near_gaussian = TRUE))
  }

  converged <- FALSE
  iter      <- 0L

  for (i in seq_len(max_iter)) {
    iter <- as.integer(i)
    eps  <- as.numeric(Y - X %*% B)

    # Adaptive mode: re-estimate kappa from current residuals
    if (adaptive) {
      eps_c <- eps - mean(eps)
      m2_c  <- mean(eps_c^2)
      m4_c  <- mean(eps_c^4)
      m6_c  <- mean(eps_c^6)
      d     <- m4_c - 3 * m2_c^2
      if (abs(d) > .Machine$double.eps * 1e6) {
        kappa <- (m6_c - 3 * m4_c * m2_c) / d
      }
    }

    # Score: Z = X'(eps * (kappa - eps^2))
    Z1 <- eps * (kappa - eps^2)
    Z  <- as.numeric(crossprod(X, Z1))

    # Jacobian: J = X' diag(3*eps^2 - kappa) X
    w  <- 3 * eps^2 - kappa
    J  <- crossprod(X * w, X)

    # Regularise if near-singular
    if (rcond(J) < 1e-12) J <- J + diag(1e-8, P)

    # Solve
    delta <- tryCatch(as.numeric(solve(J, Z)),
                      error = function(e) rep(NA_real_, P))
    if (anyNA(delta)) break

    # Step-size limiting
    step_norm <- sqrt(sum(delta^2))
    if (step_norm > step_max) delta <- delta * step_max / step_norm

    B <- B - delta

    if (step_norm < tol) {
      converged <- TRUE
      break
    }
  }

  # Divergence guard: revert to OLS if too far
  if (sqrt(sum((B - B_ols)^2)) > 10 || anyNA(B)) {
    B         <- B_ols
    converged <- FALSE
  }

  list(B = B, converged = converged, iter = iter,
       kappa = kappa, near_gaussian = FALSE)
}


#' PMM3 Newton-Raphson solver for nonlinear models
#'
#' Uses a residual function and numerical Jacobian (via \code{numDeriv})
#' instead of a fixed design matrix. This is necessary for ARIMA and other
#' models where the linearized design matrix does not capture full dynamics.
#'
#' The PMM3 estimating equation is:
#' \deqn{\sum_t \varepsilon_t(\kappa - \varepsilon_t^2) \cdot
#'       \frac{\partial \varepsilon_t}{\partial \theta} = 0}
#'
#' @param theta_init Initial parameter estimates (from MLE/CSS)
#' @param fn_residuals Function(theta) returning residual vector
#' @param kappa Moment ratio (scalar)
#' @param adaptive Re-estimate kappa each iteration? (default FALSE)
#' @param tol Convergence tolerance (default 1e-6)
#' @param max_iter Maximum NR iterations (default 100)
#' @param step_max Maximum step size (default 5.0)
#'
#' @return List: B (estimates), converged, iter, kappa, near_gaussian
#' @keywords internal
.pmm3_nr_solver_nonlinear <- function(theta_init, fn_residuals, kappa,
                                       adaptive = FALSE,
                                       tol      = 1e-6,
                                       max_iter = 100,
                                       step_max = 5.0) {

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for ARIMA PMM3. Please install it.",
         call. = FALSE)
  }

  P <- length(theta_init)
  B <- theta_init

  # Near-Gaussian: kappa undefined -> return initial
  if (is.na(kappa)) {
    return(list(B = theta_init, converged = TRUE, iter = 0L,
                kappa = NA_real_, near_gaussian = TRUE))
  }

  converged <- FALSE
  iter      <- 0L

  # Determine valid residual indices from initial residuals
  eps_full <- fn_residuals(B)
  valid_idx <- which(is.finite(eps_full))
  n_valid   <- length(valid_idx)

  # Wrapper that always returns fixed-length vector (only valid indices)
  fn_res_valid <- function(th) {
    r <- fn_residuals(th)
    out <- r[valid_idx]
    out[!is.finite(out)] <- 0  # safety: replace any new NAs with 0
    out
  }

  for (i in seq_len(max_iter)) {
    iter <- as.integer(i)

    # Current residuals (fixed length)
    eps <- fn_res_valid(B)

    # Adaptive mode: re-estimate kappa
    if (adaptive) {
      eps_c <- eps - mean(eps)
      m2_c  <- mean(eps_c^2)
      m4_c  <- mean(eps_c^4)
      m6_c  <- mean(eps_c^6)
      d     <- m4_c - 3 * m2_c^2
      if (abs(d) > .Machine$double.eps * 1e6) {
        kappa <- (m6_c - 3 * m4_c * m2_c) / d
      }
    }

    # Numerical Jacobian: J_res[i,j] = d(eps_i)/d(theta_j)
    J_res <- numDeriv::jacobian(fn_res_valid, B, method = "Richardson")

    # Score: gradient = J_res' * eps*(kappa - eps^2)
    Z1       <- eps * (kappa - eps^2)
    gradient <- as.numeric(crossprod(J_res, Z1))

    # Hessian approximation: H = J_res' * diag(3*eps^2 - kappa) * J_res
    w <- 3 * eps^2 - kappa
    H <- crossprod(J_res * w, J_res)

    # Regularise
    if (rcond(H) < 1e-12) H <- H + diag(1e-8, P)

    # Solve for delta
    delta <- tryCatch(as.numeric(solve(H, gradient)),
                      error = function(e) rep(NA_real_, P))
    if (anyNA(delta)) break

    # Step-size limiting
    step_norm <- sqrt(sum(delta^2))
    if (step_norm > step_max) delta <- delta * step_max / step_norm

    # NB: + not - because J_res = d(eps)/d(theta), and the design matrix
    # analog D = -J_res. The NR step is theta - H^{-1}*D'*Z1,
    # which equals theta + H^{-1}*J_res'*Z1 = theta + delta.
    B <- B + delta

    if (step_norm < tol) {
      converged <- TRUE
      break
    }
  }

  # Divergence guard
  if (sqrt(sum((B - theta_init)^2)) > 10 || anyNA(B)) {
    B         <- theta_init
    converged <- FALSE
  }

  list(B = B, converged = converged, iter = iter,
       kappa = kappa, near_gaussian = FALSE)
}
