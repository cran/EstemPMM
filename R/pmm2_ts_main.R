# pmm2_ts_main.R - Unified module for time series models with PMM2

#' Fit a time series model using the PMM2 method
#'
#' @param x Numeric vector of time series data
#' @param order Model order specification:
#'        - For AR models: a single integer (AR order)
#'        - For MA models: a single integer (MA order)
#'        - For ARMA models: vector c(p, q) (AR and MA orders)
#'        - For ARIMA models: vector c(p, d, q) (AR, differencing, and MA orders)
#' @param model_type String specifying the model type: "ar", "ma", "arma", or "arima"
#' @param method String: estimation method, one of "pmm2" (default), "css", "ml", "yw", "ols"
#' @param max_iter Integer: maximum number of iterations for the algorithm
#' @param tol Numeric: tolerance for convergence
#' @param include.mean Logical: whether to include a mean (intercept) term
#' @param initial List or vector of initial parameter estimates (optional)
#' @param na.action Function for handling missing values, default is na.fail
#' @param regularize Logical, add small values to diagonal for numerical stability
#' @param reg_lambda Regularization parameter (if regularize=TRUE)
#' @param verbose Logical: whether to print progress information
#'
#' @details
#' The PMM2 algorithm works as follows:
#'
#' 1. Fits an initial model using a standard method (OLS, Yule-Walker, CSS or ML)
#' 2. Computes central moments (m2, m3, m4) from initial residuals/innovations
#' 3. Uses these moments with a specialized solver (pmm2_algorithm) to find
#'    robust parameter estimates
#'
#' @return An S4 object \code{TS2fit} of the corresponding subclass
#' @export
ts_pmm2 <- function(x, order,
                    model_type = c("ar", "ma", "arma", "arima"),
                    method      = "pmm2",
                    max_iter    = 50,
                    tol         = 1e-6,
                    include.mean= TRUE,
                    initial     = NULL,
                    na.action   = na.fail,
                    regularize  = TRUE,
                    reg_lambda  = 1e-8,
                    verbose     = FALSE) {

  model_type <- match.arg(model_type)
  cl <- match.call()

  if (!is.null(na.action)) {
    x <- na.action(x)
  }

  # 1) Validate input data
  model_params <- validate_ts_parameters(x, order, model_type, include.mean)

  if (model_params$model_type == "ma") {
    q <- model_params$ma_order
    css_fit <- ma_css_fit(model_params$original_x, q, include.mean, verbose)

    if (method == "css") {
      moments <- compute_moments(css_fit$residuals)
      return(new("MAPMM2",
                 coefficients    = as.numeric(css_fit$coefficients),
                 residuals       = as.numeric(css_fit$residuals),
                 m2              = as.numeric(moments$m2),
                 m3              = as.numeric(moments$m3),
                 m4              = as.numeric(moments$m4),
                 convergence     = css_fit$convergence,
                 iterations      = as.numeric(css_fit$iterations),
                 call            = cl,
                 model_type      = "ma",
                 intercept       = if (include.mean) as.numeric(css_fit$intercept) else 0,
                 original_series = as.numeric(model_params$original_x),
                 order           = list(ar = 0L, ma = q, d = 0L)))
    }

    pmm2_fit <- ma_pmm2_fit(model_params$original_x, q, css_fit,
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)
    moments <- compute_moments(pmm2_fit$innovations)

    return(new("MAPMM2",
               coefficients    = as.numeric(pmm2_fit$coefficients),
               residuals       = as.numeric(pmm2_fit$innovations),
               m2              = as.numeric(moments$m2),
               m3              = as.numeric(moments$m3),
               m4              = as.numeric(moments$m4),
               convergence     = pmm2_fit$convergence,
               iterations      = as.numeric(pmm2_fit$iterations),
               call            = cl,
               model_type      = "ma",
               intercept       = if (include.mean) as.numeric(pmm2_fit$intercept) else 0,
               original_series = as.numeric(model_params$original_x),
               order           = list(ar = 0L, ma = q, d = 0L)))
  }

  if (model_params$model_type == "arima") {
    p <- model_params$ar_order
    d <- model_params$d
    q <- model_params$ma_order

    css_fit <- tryCatch(
      stats::arima(model_params$original_x,
                   order = c(p, d, q),
                   method = "CSS-ML",
                   include.mean = include.mean),
      error = function(e) NULL
    )

    if (is.null(css_fit)) {
      stop("Failed to estimate ARIMA model using classical method")
    }

    coef_names <- names(css_fit$coef)
    ar_css <- if (p > 0) as.numeric(css_fit$coef[paste0("ar", seq_len(p))]) else numeric(0)
    ma_css <- if (q > 0) as.numeric(css_fit$coef[paste0("ma", seq_len(q))]) else numeric(0)

    intercept_css <- 0
    if (include.mean && d == 0) {
      intercept_name <- setdiff(coef_names,
                                c(paste0("ar", seq_len(p)), paste0("ma", seq_len(q))))
      if (length(intercept_name) > 0) {
        intercept_css <- as.numeric(css_fit$coef[intercept_name[1]])
      }
    }

    residuals_css <- as.numeric(css_fit$residuals)
    residuals_css[is.na(residuals_css)] <- 0

    x_diff <- if (d > 0) diff(model_params$original_x, differences = d)
              else model_params$original_x
    res_diff <- tail(residuals_css, length(x_diff))

    include_intercept_diff <- include.mean && d == 0

    if (method == "css") {
      res_clean <- res_diff[is.finite(res_diff)]
      moments <- compute_moments(res_clean)
      return(new("ARIMAPMM2",
                 coefficients    = as.numeric(c(ar_css, ma_css)),
                 residuals       = as.numeric(residuals_css),
                 m2              = as.numeric(moments$m2),
                 m3              = as.numeric(moments$m3),
                 m4              = as.numeric(moments$m4),
                 convergence     = TRUE,
                 iterations      = 1L,
                 call            = cl,
                 model_type      = "arima",
                 intercept       = if (include_intercept_diff) intercept_css else 0,
                 original_series = as.numeric(model_params$original_x),
                 order           = list(ar = p, ma = q, d = d)))
    }

    design <- arma_build_design(x_diff, res_diff,
                                 p = p, q = q,
                                 intercept = intercept_css,
                                 include_intercept = include_intercept_diff)

    moments <- compute_moments(res_diff[is.finite(res_diff)])
    b_init <- c(if (include_intercept_diff) 0 else NULL, ar_css, ma_css)

    algo_res <- pmm2_algorithm(
      b_init = b_init,
      X = design$X,
      y = design$y,
      m2 = moments$m2,
      m3 = moments$m3,
      m4 = moments$m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    if (include_intercept_diff) {
      intercept_hat <- algo_res$b[1]
      ar_hat <- if (p > 0) algo_res$b[1 + seq_len(p)] else numeric(0)
      ma_hat <- if (q > 0) algo_res$b[1 + p + seq_len(q)] else numeric(0)
    } else {
      intercept_hat <- 0
      ar_hat <- if (p > 0) algo_res$b[seq_len(p)] else numeric(0)
      ma_hat <- if (q > 0) algo_res$b[p + seq_len(q)] else numeric(0)
    }

    x_mean <- if (include_intercept_diff) intercept_hat else 0
    x_centered <- x_diff - x_mean

    model_info <- list(
      ar_order = p,
      ma_order = q,
      d = d,
      model_type = "arima",
      include.mean = include.mean,
      innovations = res_diff,
      x = x_centered,
      x_mean = x_mean,
      original_x = model_params$original_x,
      verbose = verbose
    )

    final_coef <- c(ar_hat, ma_hat)
    final_res <- compute_ts_residuals(final_coef, model_info)
    res_clean <- final_res[is.finite(final_res)]
    moments_final <- compute_moments(res_clean)

    return(new("ARIMAPMM2",
               coefficients    = as.numeric(final_coef),
               residuals       = as.numeric(final_res),
               m2              = as.numeric(moments_final$m2),
               m3              = as.numeric(moments_final$m3),
               m4              = as.numeric(moments_final$m4),
               convergence     = algo_res$convergence,
               iterations      = as.numeric(algo_res$iterations),
               call            = cl,
               model_type      = "arima",
               intercept       = if (include_intercept_diff) as.numeric(intercept_hat) else 0,
               original_series = as.numeric(model_params$original_x),
               order           = list(ar = p, ma = q, d = d)))
  }

  # 2) Obtain initial estimates
  init <- get_initial_estimates(model_params, initial, method, verbose)
  b_init      <- init$b_init
  x_mean      <- init$x_mean
  innovations <- init$innovations
  x_centered  <- init$x_centered
  orig_x      <- init$orig_x
  m2          <- init$m2
  m3          <- init$m3
  m4          <- init$m4

  # 3) Create design matrix
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

  # 4) If method == "pmm2", use the PMM2 algorithm
  if (method == "pmm2") {
    if (verbose) cat("Starting PMM2 optimization...\n")

    model_info <- list(
      ar_order = model_params$ar_order,
      ma_order = model_params$ma_order,
      d = model_params$d,
      model_type = model_params$model_type,
      include.mean = model_params$include.mean,
      innovations = innovations,
      x = x_centered,
      x_mean = x_mean,
      original_x = orig_x,
      verbose = verbose
    )

    result <- pmm2_algorithm(
      b_init = b_init,
      X = dm$X,
      y = dm$y,
      m2 = m2,
      m3 = m3,
      m4 = m4,
      max_iter = max_iter,
      tol = tol,
      regularize = regularize,
      reg_lambda = reg_lambda,
      verbose = verbose
    )

    final_coef <- result$b
    converged <- result$convergence
    iterations <- result$iterations

    # Compute final residuals
    final_res <- compute_ts_residuals(final_coef, model_info)
  } else {
    # For other methods, simply use initial estimates
    final_coef <- b_init
    converged <- TRUE
    iterations <- 0
    final_res <- innovations
  }

  # 5) Create the appropriate class object
  if (model_type == "ar") {
    result_class <- "ARPMM2"
  } else if (model_type == "ma") {
    result_class <- "MAPMM2"
  } else if (model_type == "arma") {
    result_class <- "ARMAPMM2"
  } else if (model_type == "arima") {
    result_class <- "ARIMAPMM2"
  } else {
    result_class <- "TS2fit"  # Base class by default
  }

  # Create and return the object of the appropriate class
  new(result_class,
      coefficients    = as.numeric(final_coef),
      residuals       = as.numeric(final_res),
      m2              = as.numeric(m2),
      m3              = as.numeric(m3),
      m4              = as.numeric(m4),
      convergence     = converged,
      iterations      = as.numeric(iterations),
      call            = cl,
      model_type      = model_type,
      intercept       = as.numeric(x_mean),
      original_series = as.numeric(orig_x),
      order           = list(ar = model_params$ar_order,
                             ma = model_params$ma_order,
                             d = model_params$d))
}

#' Fit an AR model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{ARPMM2} containing fitted autoregressive
#'   coefficients, residuals, central moment estimates (m2-m4), model order,
#'   intercept, original series, and convergence diagnostics.
#' @export
ar_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ar", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Fit an MA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{MAPMM2} containing moving-average
#'   coefficients, residual innovations, central moments, model order,
#'   intercept, original series, and convergence diagnostics.
#' @export
ma_pmm2 <- function(x, order = 1, method = "pmm2", max_iter = 50, tol = 1e-6,
                    include.mean = TRUE, initial = NULL, na.action = na.fail,
                    regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "ma", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Fit an ARMA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{ARMAPMM2} containing fitted AR and MA
#'   coefficients, residuals, central moments, model specification, intercept,
#'   original series, and convergence diagnostics.
#' @export
arma_pmm2 <- function(x, order = c(1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                      include.mean = TRUE, initial = NULL, na.action = na.fail,
                      regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arma", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

#' Fit an ARIMA model using PMM2 (wrapper)
#'
#' @inheritParams ts_pmm2
#' @return An S4 object of class \code{ARIMAPMM2} containing fitted AR and MA
#'   coefficients, residual series, central moments, differencing order,
#'   intercept, original series, and convergence diagnostics.
#' @export
arima_pmm2 <- function(x, order = c(1, 1, 1), method = "pmm2", max_iter = 50, tol = 1e-6,
                       include.mean = TRUE, initial = NULL, na.action = na.fail,
                       regularize = TRUE, reg_lambda = 1e-8, verbose = FALSE) {
  ts_pmm2(x, order = order, model_type = "arima", method = method,
          max_iter = max_iter, tol = tol,
          include.mean = include.mean, initial = initial,
          na.action = na.action, regularize = regularize,
          reg_lambda = reg_lambda, verbose = verbose)
}

# --- MA utilities ---------------------------------------------------------

ma_css_fit <- function(x, q, include_mean = TRUE, verbose = FALSE) {
  fit <- tryCatch(
    stats::arima(x, order = c(0, 0, q), method = "CSS-ML",
                 include.mean = include_mean),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    if (verbose) {
      cat("Failed to estimate MA model via stats::arima: returning zero coefficients\n")
    }
    coef_css <- rep(0, q)
    intercept <- if (include_mean) mean(x) else 0
    residuals <- x - intercept
    residuals[is.na(residuals)] <- 0
    return(list(
      coefficients = coef_css,
      intercept = intercept,
      residuals = residuals,
      convergence = FALSE,
      iterations = 0L
    ))
  }

  names_coef <- names(fit$coef)
  coef_css <- numeric(q)
  for (j in seq_len(q)) {
    coef_name <- paste0("ma", j)
    if (coef_name %in% names_coef) {
      coef_css[j] <- fit$coef[coef_name]
    } else if (length(fit$coef) >= j) {
      coef_css[j] <- fit$coef[j]
    } else {
      coef_css[j] <- 0
    }
  }

  intercept <- 0
  if (include_mean) {
    intercept_name <- setdiff(names_coef, paste0("ma", seq_len(q)))
    if (length(intercept_name) > 0) {
      intercept <- as.numeric(fit$coef[intercept_name[1]])
    }
  }

  residuals <- as.numeric(fit$residuals)
  residuals[is.na(residuals)] <- 0

  list(
    coefficients = coef_css,
    intercept = intercept,
    residuals = residuals,
    convergence = TRUE,
    iterations = 1L
  )
}

ma_pmm2_fit <- function(x, q, css_fit, max_iter = 50, tol = 1e-6, verbose = FALSE) {
  design <- ma_build_design(css_fit$intercept, css_fit$residuals, x, q)
  moments <- compute_moments(css_fit$residuals)

  b_init <- c(0, css_fit$coefficients)
  solve_res <- ma_solve_pmm2(b_init, design$X, design$y,
                             moments$m2, moments$m3, moments$m4,
                             max_iter = max_iter, tol = tol,
                             verbose = verbose)

  theta <- solve_res$coefficients
  innovations <- ma_compute_innovations(x - css_fit$intercept, theta, q)

  list(
    coefficients = theta,
    intercept = css_fit$intercept,
    innovations = innovations,
    convergence = solve_res$convergence,
    iterations = solve_res$iterations
  )
}

ma_build_design <- function(intercept, residuals, x, q) {
  idx <- seq.int(q + 1L, length(x))
  X <- matrix(1, nrow = length(idx), ncol = q + 1L)
  for (j in seq_len(q)) {
    X[, j + 1L] <- residuals[idx - j]
  }
  y <- x[idx] - intercept
  list(X = X, y = y)
}

ma_solve_pmm2 <- function(b_init, X, Y, m2, m3, m4,
                          max_iter = 50, tol = 1e-6,
                          verbose = FALSE) {
  b <- as.numeric(b_init)
  iterations <- 0L
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    iterations <- iter
    S <- as.vector(X %*% b)
    Z1 <- m3 * S^2 + (m4 - m2^2 - 2 * m3 * Y) * S +
      (m3 * Y^2 - (m4 - m2^2) * Y - m2 * m3)
    Z <- as.numeric(t(X) %*% Z1)
    JZ11 <- 2 * m3 * S + (m4 - m2^2 - 2 * m3 * Y)
    J <- t(X) %*% (X * JZ11)
    step <- tryCatch(solve(J, Z), error = function(e) NULL)
    if (is.null(step)) {
      if (verbose) cat("System is singular at iteration", iter, "\n")
      break
    }
    b_new <- b - step
    if (sqrt(sum((b_new - b)^2)) < tol) {
      b <- b_new
      converged <- TRUE
      break
    }
    b <- b_new
  }
  list(coefficients = b[-1], convergence = converged, iterations = iterations)
}

ma_compute_innovations <- function(x, theta, q) {
  n <- length(x)
  innovations <- numeric(n)
  history <- rep(0, q)
  for (t in seq_len(n)) {
    ma_component <- if (q > 0) sum(theta * history) else 0
    innovations[t] <- x[t] - ma_component
    if (q > 0) {
      history <- c(innovations[t], history)[seq_len(q)]
    }
  }
  innovations
}

arma_build_design <- function(x, residuals, p, q, intercept = 0, include_intercept = FALSE) {
  n <- length(x)
  max_lag <- max(p, q)
  if (n <= max_lag) {
    stop("Insufficient data to build ARMA design matrix")
  }

  idx <- seq.int(max_lag + 1L, n)
  columns <- list()
  if (include_intercept) {
    columns <- c(columns, list(rep(1, length(idx))))
  }
  if (p > 0) {
    for (j in seq_len(p)) {
      columns <- c(columns, list(x[idx - j]))
    }
  }
  if (q > 0) {
    for (j in seq_len(q)) {
      columns <- c(columns, list(residuals[idx - j]))
    }
  }

  X <- if (length(columns) > 0) do.call(cbind, columns) else matrix(0, length(idx), 0)
  y <- x[idx] - intercept
  list(X = X, y = y)
}
