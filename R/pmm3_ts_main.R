# pmm3_ts_main.R - Time series models with PMM3 (S=3)

#' Fit a time series model using PMM3
#'
#' Core function that fits AR, MA, ARMA, or ARIMA models using the
#' Polynomial Maximization Method of order 3 (PMM3). Designed for
#' symmetric platykurtic innovations.
#'
#' @param x Numeric vector of time series data
#' @param order Model order specification (see \code{\link{ts_pmm2}} for format)
#' @param model_type Character: "ar", "ma", "arma", or "arima"
#' @param max_iter Integer: maximum NR iterations (default 100)
#' @param tol Numeric: convergence tolerance (default 1e-6)
#' @param adaptive Logical: re-estimate kappa each iteration (default FALSE)
#' @param step_max Numeric: maximum NR step size (default 5.0)
#' @param include.mean Logical: include mean/intercept term (default TRUE)
#' @param initial Optional initial parameter estimates
#' @param na.action Function for handling missing values (default na.fail)
#' @param verbose Logical: print progress information (default FALSE)
#'
#' @details
#' The PMM3 time series algorithm:
#' 1. Obtains initial estimates via classical methods (OLS/YW for AR, CSS for MA/ARMA/ARIMA)
#' 2. Computes moments m2, m4, m6 and kappa from initial residuals
#' 3. Checks symmetry: warns if |gamma3| > 0.3
#' 4. Applies Newton-Raphson with PMM3 score equations
#'
#' PMM3 is beneficial when innovations are symmetric and platykurtic
#' (gamma4 < 0), e.g. uniform, truncated normal.
#'
#' @return An S4 object of the appropriate TS3fit subclass
#' @export
ts_pmm3 <- function(x, order,
                    model_type = c("ar", "ma", "arma", "arima"),
                    max_iter = 100, tol = 1e-6,
                    adaptive = FALSE, step_max = 5.0,
                    include.mean = TRUE, initial = NULL,
                    na.action = na.fail, verbose = FALSE) {

  model_type <- match.arg(model_type)
  cl <- match.call()

  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # 1) Validate input
  model_params <- validate_ts_parameters(x, order, model_type, include.mean)

  # 2) Get initial estimates (reuse PMM2 infrastructure)
  init <- get_initial_estimates(model_params, initial, method = "pmm2", verbose = verbose)
  b_init     <- init$b_init
  x_mean     <- init$x_mean
  innovations <- init$innovations
  x_centered <- init$x_centered
  orig_x     <- init$orig_x

  # 3) Compute PMM3 moments from initial residuals
  mom <- compute_moments_pmm3(innovations)

  if (verbose) {
    cat("PMM3 moments from initial residuals:\n")
    cat("  m2 =", mom$m2, "  m4 =", mom$m4, "  m6 =", mom$m6, "\n")
    cat("  gamma3 =", mom$gamma3, "  gamma4 =", mom$gamma4, "\n")
    cat("  kappa =", mom$kappa, "\n")
  }

  # Symmetry warning
  if (abs(mom$gamma3) > 0.3) {
    warning(sprintf(
      "Innovations appear asymmetric (|gamma3| = %.3f > 0.3). Consider using PMM2 time series functions instead.",
      abs(mom$gamma3)
    ))
  }

  # 4-5) Apply PMM3 solver (method depends on model type)
  if (verbose) cat("Starting PMM3 Newton-Raphson...\n")

  if (model_type == "arima" && model_params$d > 0) {
    # ARIMA: use nonlinear solver with stats::arima for proper residual computation
    arima_order <- c(model_params$ar_order, model_params$d, model_params$ma_order)
    inc_mean_arima <- model_params$include.mean && model_params$d == 0

    fn_residuals <- function(theta) {
      ar_p <- model_params$ar_order
      ma_q <- model_params$ma_order
      fixed_params <- theta
      if (inc_mean_arima) {
        fixed_params <- c(fixed_params, x_mean)
      }
      fit_tmp <- tryCatch(
        stats::arima(orig_x, order = arima_order, fixed = fixed_params,
                     include.mean = inc_mean_arima, transform.pars = FALSE),
        error = function(e) NULL
      )
      if (is.null(fit_tmp)) return(rep(NA_real_, length(orig_x)))
      as.numeric(fit_tmp$residuals)
    }

    nr <- .pmm3_nr_solver_nonlinear(
      theta_init = b_init,
      fn_residuals = fn_residuals,
      kappa    = mom$kappa,
      adaptive = adaptive,
      tol      = tol,
      max_iter = max_iter,
      step_max = step_max
    )
  } else {
    # AR/MA/ARMA: use linearized design matrix solver (fast, works well)
    dm <- create_ts_design_matrix(
      x = x_centered,
      model_info = list(
        ar_order = model_params$ar_order,
        ma_order = model_params$ma_order,
        d = model_params$d,
        model_type = model_params$model_type,
        include.mean = model_params$include.mean
      ),
      innovations = innovations
    )

    nr <- .pmm3_nr_solver(
      B_ols    = b_init,
      X        = dm$X,
      Y        = dm$y,
      kappa    = mom$kappa,
      adaptive = adaptive,
      tol      = tol,
      max_iter = max_iter,
      step_max = step_max
    )
  }

  final_coef <- nr$B
  converged  <- nr$converged
  iterations <- nr$iter

  # 6) Compute final residuals
  model_info <- list(
    ar_order   = model_params$ar_order,
    ma_order   = model_params$ma_order,
    d          = model_params$d,
    model_type = model_params$model_type,
    include.mean = model_params$include.mean,
    innovations = innovations,
    x          = x_centered,
    x_mean     = x_mean,
    original_x = orig_x,
    verbose    = verbose
  )

  final_res <- compute_ts_residuals(final_coef, model_info)
  res_clean <- final_res[is.finite(final_res)]
  mom_final <- compute_moments_pmm3(res_clean)
  vf <- pmm3_variance_factor(mom$m2, mom$m4, mom$m6)

  if (verbose) {
    cat("PMM3 completed. Converged:", converged, " Iterations:", iterations, "\n")
  }

  # 7) Determine result class
  result_class <- switch(model_type,
    ar    = "ARPMM3",
    ma    = "MAPMM3",
    arma  = "ARMAPMM3",
    arima = "ARIMAPMM3",
    "TS3fit"
  )

  new(result_class,
    coefficients    = as.numeric(final_coef),
    residuals       = as.numeric(final_res),
    m2              = as.numeric(mom$m2),
    m4              = as.numeric(mom$m4),
    m6              = as.numeric(mom$m6),
    gamma4          = as.numeric(vf$gamma4),
    gamma6          = as.numeric(vf$gamma6),
    g_coefficient   = as.numeric(if (is.na(vf$g3)) 1.0 else vf$g3),
    kappa           = as.numeric(if (is.na(mom$kappa)) NA_real_ else mom$kappa),
    convergence     = converged,
    iterations      = as.numeric(iterations),
    call            = cl,
    model_type      = model_type,
    intercept       = as.numeric(x_mean),
    original_series = as.numeric(orig_x),
    order           = list(
      ar = model_params$ar_order,
      ma = model_params$ma_order,
      d  = model_params$d
    )
  )
}


# =============================================================================
# Wrapper functions
# =============================================================================

#' Fit an AR model using PMM3
#'
#' Estimates autoregressive model parameters using the Polynomial Maximization
#' Method of order 3 (PMM3). Designed for symmetric platykurtic innovations.
#'
#' @param x Numeric vector of time series data
#' @param order Integer: AR order (default 1)
#' @param max_iter Integer: maximum NR iterations (default 100)
#' @param tol Numeric: convergence tolerance (default 1e-6)
#' @param adaptive Logical: re-estimate kappa each iteration (default FALSE)
#' @param step_max Numeric: maximum NR step size (default 5.0)
#' @param include.mean Logical: include mean term (default TRUE)
#' @param initial Optional initial parameter estimates
#' @param na.action Function for handling missing values (default na.fail)
#' @param verbose Logical: print progress information (default FALSE)
#'
#' @return An S4 object of class \code{ARPMM3}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- arima.sim(n = 200, list(ar = 0.7),
#'                innov = runif(200, -sqrt(3), sqrt(3)))
#' fit <- ar_pmm3(x, order = 1)
#' coef(fit)
#' }
#'
#' @seealso \code{\link{ar_pmm2}}, \code{\link{ma_pmm3}}, \code{\link{lm_pmm3}}
#' @export
ar_pmm3 <- function(x, order = 1, max_iter = 100, tol = 1e-6,
                    adaptive = FALSE, step_max = 5.0,
                    include.mean = TRUE, initial = NULL,
                    na.action = na.fail, verbose = FALSE) {
  ts_pmm3(x, order = order, model_type = "ar",
          max_iter = max_iter, tol = tol,
          adaptive = adaptive, step_max = step_max,
          include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}


#' Fit an MA model using PMM3
#'
#' Estimates moving-average model parameters using PMM3.
#'
#' @inheritParams ar_pmm3
#'
#' @return An S4 object of class \code{MAPMM3}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- arima.sim(n = 200, list(ma = 0.6),
#'                innov = runif(200, -sqrt(3), sqrt(3)))
#' fit <- ma_pmm3(x, order = 1)
#' coef(fit)
#' }
#'
#' @seealso \code{\link{ma_pmm2}}, \code{\link{ar_pmm3}}, \code{\link{lm_pmm3}}
#' @export
ma_pmm3 <- function(x, order = 1, max_iter = 100, tol = 1e-6,
                    adaptive = FALSE, step_max = 5.0,
                    include.mean = TRUE, initial = NULL,
                    na.action = na.fail, verbose = FALSE) {
  ts_pmm3(x, order = order, model_type = "ma",
          max_iter = max_iter, tol = tol,
          adaptive = adaptive, step_max = step_max,
          include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}


#' Fit an ARMA model using PMM3
#'
#' Estimates ARMA model parameters using PMM3.
#'
#' @inheritParams ar_pmm3
#' @param order Numeric vector c(p, q): AR and MA orders
#'
#' @return An S4 object of class \code{ARMAPMM3}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- arima.sim(n = 250, list(ar = 0.7, ma = -0.3),
#'                innov = runif(250, -sqrt(3), sqrt(3)))
#' fit <- arma_pmm3(x, order = c(1, 1))
#' coef(fit)
#' }
#'
#' @seealso \code{\link{arma_pmm2}}, \code{\link{arima_pmm3}}, \code{\link{lm_pmm3}}
#' @export
arma_pmm3 <- function(x, order = c(1, 1), max_iter = 100, tol = 1e-6,
                      adaptive = FALSE, step_max = 5.0,
                      include.mean = TRUE, initial = NULL,
                      na.action = na.fail, verbose = FALSE) {
  ts_pmm3(x, order = order, model_type = "arma",
          max_iter = max_iter, tol = tol,
          adaptive = adaptive, step_max = step_max,
          include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}


#' Fit an ARIMA model using PMM3
#'
#' Estimates ARIMA model parameters using PMM3.
#'
#' @inheritParams ar_pmm3
#' @param order Numeric vector c(p, d, q): AR, differencing, and MA orders
#'
#' @return An S4 object of class \code{ARIMAPMM3}
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(arima.sim(n = 200, list(ar = 0.6),
#'             innov = runif(200, -sqrt(3), sqrt(3))))
#' fit <- arima_pmm3(x, order = c(1, 1, 0))
#' coef(fit)
#' }
#'
#' @seealso \code{\link{arima_pmm2}}, \code{\link{arma_pmm3}}, \code{\link{lm_pmm3}}
#' @export
arima_pmm3 <- function(x, order = c(1, 1, 1), max_iter = 100, tol = 1e-6,
                       adaptive = FALSE, step_max = 5.0,
                       include.mean = TRUE, initial = NULL,
                       na.action = na.fail, verbose = FALSE) {
  ts_pmm3(x, order = order, model_type = "arima",
          max_iter = max_iter, tol = tol,
          adaptive = adaptive, step_max = step_max,
          include.mean = include.mean, initial = initial,
          na.action = na.action, verbose = verbose)
}
