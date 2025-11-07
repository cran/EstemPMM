# pmm2_monte_carlo.R - Monte Carlo comparison of PMM2 method accuracy

.pmm2_default_burn_in <- 100L

#' Monte Carlo comparison of PMM2 estimation methods
#'
#' Function generates time series for given models, repeatedly estimates
#' parameters using different methods and compares their accuracy by MSE criterion.
#' Additionally outputs theoretical and empirical characteristics of the innovation
#' distribution (skewness, excess kurtosis, theoretical gain of PMM2).
#'
#' @param model_specs List of model specifications. Each element must contain:
#'   \describe{
#'     \item{model}{"ar", "ma" or "arma"}
#'     \item{order}{order (for AR/MA) or vector c(p, q) for ARMA}
#'     \item{theta}{numeric vector of true parameters; for ARMA a list
#'                  `list(ar = ..., ma = ...)`}
#'     \item{label}{(optional) model name in report}
#'     \item{innovations}{(optional) description of innovation distribution:
#'           \code{list(type = "gamma", shape = 2)},
#'           \code{list(type = "student_t", df = 5)}, etc. Can also
#'           pass an arbitrary generation function via \code{generator}.}
#'   }
#' @param methods Vector of estimation methods (e.g., `c("css","pmm2")`).
#'                The first method is considered baseline for relative MSE calculation.
#' @param n Sample size for simulation.
#' @param n_sim Number of Monte Carlo experiments.
#' @param innovations Function or distribution description, used by
#'                    default for all models (if not specified in spec).
#' @param seed Initial seed for random number generator (optional).
#' @param include.mean Logical flag: whether to include intercept during estimation.
#' @param progress Logical flag: print Monte Carlo progress.
#' @param verbose Whether to print diagnostic messages on failures.
#'
#' @return List with three components:
#'   \describe{
#'     \item{parameter_results}{MSE and relative MSE for each parameter}
#'     \item{summary}{Averaged MSE over parameters for each model/method}
#'     \item{gain}{Comparison of theoretical and empirical PMM2 gain}
#'   }
#'
#' @export
pmm2_monte_carlo_compare <- function(model_specs,
                                     methods = c("css", "pmm2"),
                                     n,
                                     n_sim,
                                     innovations = list(type = "gaussian"),
                                     seed = NULL,
                                     include.mean = TRUE,
                                     progress = interactive(),
                                     verbose = FALSE) {
  if (missing(model_specs) || length(model_specs) == 0L) {
    stop("model_specs must contain at least one model")
  }
  if (missing(n) || n <= 0) {
    stop("n must be positive")
  }
  if (missing(n_sim) || n_sim <= 0) {
    stop("n_sim must be positive")
  }
  if (length(methods) < 2L) {
    warning("Less than two methods provided; relative MSE is calculated relative to the first")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  spec_results <- lapply(model_specs, function(spec) {
    innov_info <- resolve_innovation_spec(spec$innovations, innovations)
    compare_single_spec(
      spec = spec,
      methods = methods,
      n = n,
      n_sim = n_sim,
      innov_info = innov_info,
      include.mean = include.mean,
      progress = progress,
      verbose = verbose
    )
  })

  parameter_results <- do.call(rbind, lapply(spec_results, `[[`, "parameter"))
  summary_results <- do.call(rbind, lapply(spec_results, `[[`, "summary"))
  gain_results <- do.call(rbind, lapply(spec_results, `[[`, "gain"))

  list(
    parameter_results = parameter_results,
    summary = summary_results,
    gain = gain_results
  )
}

compare_single_spec <- function(spec,
                                methods,
                                n,
                                n_sim,
                                innov_info,
                                include.mean,
                                progress,
                                verbose) {
  model_type <- match.arg(spec$model, c("ar", "ma", "arma", "arima"))
  label <- if (!is.null(spec$label)) spec$label else default_label(spec)

  true_theta <- extract_true_params(spec)
  n_params <- length(true_theta)

  sum_sq <- matrix(0, nrow = length(methods), ncol = n_params)
  success <- integer(length(methods))
  collected_innov <- numeric(0)

  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)
  }

  for (sim in seq_len(n_sim)) {
    series_obj <- tryCatch(
      simulate_series(spec, n, innov_info),
      error = function(e) {
        if (verbose) {
          message("Series generation error: ", conditionMessage(e))
        }
        NULL
      }
    )

    if (is.null(series_obj)) {
      if (progress) utils::setTxtProgressBar(pb, sim)
      next
    }

    collected_innov <- c(collected_innov, series_obj$innovations)

    for (m_idx in seq_along(methods)) {
      method <- methods[m_idx]
      est <- tryCatch(
        fit_and_extract(series_obj$series, spec, method, include.mean),
        error = function(e) {
          if (verbose) {
            message("Estimation error (", label, ", method=", method,
                    "): ", conditionMessage(e))
          }
          rep(NA_real_, n_params)
        }
      )

      if (!anyNA(est) && length(est) == n_params) {
        diff_sq <- (est - true_theta)^2
        sum_sq[m_idx, ] <- sum_sq[m_idx, ] + diff_sq
        success[m_idx] <- success[m_idx] + 1L
      }
    }

    if (progress) utils::setTxtProgressBar(pb, sim)
  }

  if (progress) close(pb)

  mse <- matrix(NA_real_, nrow = length(methods), ncol = n_params)
  for (m_idx in seq_along(methods)) {
    if (success[m_idx] > 0L) {
      mse[m_idx, ] <- sum_sq[m_idx, ] / success[m_idx]
    }
  }

  baseline_mse <- mse[1, ]

  empirical_moments <- if (length(collected_innov)) compute_moments(collected_innov) else list(c3 = NA_real_, c4 = NA_real_, g = NA_real_)
  theoretical_moments <- innov_info$theoretical

  parameter_names <- param_names(spec)
  mse_vec <- as.vector(t(mse))
  ratios <- as.vector(t(mse / matrix(baseline_mse, nrow = length(methods), ncol = n_params, byrow = TRUE)))

  parameter_df <- data.frame(
    model = label,
    method = rep(methods, each = n_params),
    parameter = rep(parameter_names, times = length(methods)),
    mse = mse_vec,
    mse_ratio = ratios,
    successful_runs = rep(success, each = n_params),
    total_runs = n_sim,
    theoretical_c3 = theoretical_moments$c3,
    theoretical_c4 = theoretical_moments$c4,
    theoretical_g = theoretical_moments$g,
    empirical_c3 = empirical_moments$c3,
    empirical_c4 = empirical_moments$c4,
    empirical_g = empirical_moments$g,
    stringsAsFactors = FALSE
  )

  summary_df <- aggregate(
    list(mse = parameter_df$mse),
    by = list(model = parameter_df$model,
              method = parameter_df$method),
    FUN = mean,
    na.rm = TRUE
  )

  summary_df$mse_ratio <- NA_real_
  for (mod in unique(summary_df$model)) {
    idx <- which(summary_df$model == mod)
    baseline_value <- summary_df$mse[idx[1]]
    summary_df$mse_ratio[idx] <- summary_df$mse[idx] / baseline_value
  }

  summary_df$theoretical_c3 <- theoretical_moments$c3
  summary_df$theoretical_c4 <- theoretical_moments$c4
  summary_df$theoretical_g <- theoretical_moments$g
  summary_df$empirical_c3 <- empirical_moments$c3
  summary_df$empirical_c4 <- empirical_moments$c4
  summary_df$empirical_g <- empirical_moments$g

  observed_ratio <- NA_real_
  if ("pmm2" %in% methods) {
    idx_pmm2 <- which(summary_df$method == "pmm2")[1]
    if (!is.na(idx_pmm2)) {
      observed_ratio <- summary_df$mse_ratio[idx_pmm2]
    }
  }

  gain_df <- data.frame(
    model = label,
    innovation = innov_info$label,
    baseline_method = methods[1],
    pmm2_method = if ("pmm2" %in% methods) "pmm2" else NA_character_,
    theoretical_c3 = theoretical_moments$c3,
    theoretical_c4 = theoretical_moments$c4,
    theoretical_g = theoretical_moments$g,
    empirical_c3 = empirical_moments$c3,
    empirical_c4 = empirical_moments$c4,
    empirical_g = empirical_moments$g,
    observed_ratio = observed_ratio,
    stringsAsFactors = FALSE
  )

  list(
    parameter = parameter_df,
    summary = summary_df,
    gain = gain_df
  )
}

simulate_series <- function(spec, n, innov_info) {
  model_type <- spec$model
  record_env <- new.env(parent = emptyenv())
  rand_gen <- function(nn, ...) {
    z <- innov_info$generator(nn)
    record_env$innov <- z
    z
  }

  series <- if (model_type == "ar") {
    stats::arima.sim(model = list(ar = spec$theta), n = n,
                     rand.gen = rand_gen)
  } else if (model_type == "ma") {
    stats::arima.sim(model = list(ma = spec$theta), n = n,
                     rand.gen = rand_gen)
  } else if (model_type == "arma") {
    theta <- spec$theta
    stats::arima.sim(model = list(ar = theta$ar, ma = theta$ma), n = n,
                     rand.gen = rand_gen)
  } else if (model_type == "arima") {
    theta <- spec$theta
    stats::arima.sim(model = list(order = spec$order,
                                  ar = theta$ar,
                                  ma = theta$ma),
                     n = n,
                     rand.gen = rand_gen)
  } else {
    stop("Unknown model type: ", model_type)
  }

  list(
    series = series,
    innovations = record_env$innov
  )
}

fit_and_extract <- function(series, spec, method, include.mean) {
  model_type <- spec$model
  if (model_type == "ar") {
    fit <- ar_pmm2(series, order = spec$order, method = method,
                   include.mean = include.mean, verbose = FALSE)
    est <- fit@coefficients[seq_len(spec$order)]
  } else if (model_type == "ma") {
    fit <- ma_pmm2(series, order = spec$order, method = method,
                   include.mean = include.mean, verbose = FALSE)
    est <- fit@coefficients[seq_len(spec$order)]
  } else if (model_type == "arma") {
    order <- spec$order
    fit <- arma_pmm2(series, order = order, method = method,
                     include.mean = include.mean, verbose = FALSE)
    p <- order[1]
    q <- order[2]
    est <- fit@coefficients[seq_len(p + q)]
  } else {
    order <- spec$order
    fit <- arima_pmm2(series, order = order, method = method,
                      include.mean = include.mean, verbose = FALSE)
    p <- order[1]
    q <- order[3]
    est <- fit@coefficients[seq_len(p + q)]
  }
  as.numeric(est)
}

extract_true_params <- function(spec) {
  if (spec$model == "arma" || spec$model == "arima") {
    theta <- spec$theta
    c(if (length(theta$ar)) theta$ar else numeric(0),
      if (length(theta$ma)) theta$ma else numeric(0))
  } else {
    spec$theta
  }
}

param_names <- function(spec) {
  if (spec$model == "ar") {
    paste0("ar", seq_len(spec$order))
  } else if (spec$model == "ma") {
    paste0("ma", seq_len(spec$order))
  } else if (spec$model == "arma") {
    p <- spec$order[1]
    q <- spec$order[2]
    c(if (p > 0) paste0("ar", seq_len(p)) else character(0),
      if (q > 0) paste0("ma", seq_len(q)) else character(0))
  } else {
    p <- spec$order[1]
    q <- spec$order[3]
    c(if (p > 0) paste0("ar", seq_len(p)) else character(0),
      if (q > 0) paste0("ma", seq_len(q)) else character(0))
  }
}

default_label <- function(spec) {
  if (spec$model == "ar") {
    paste0("AR(", spec$order, ")")
  } else if (spec$model == "ma") {
    paste0("MA(", spec$order, ")")
  } else if (spec$model == "arma") {
    paste0("ARMA(", spec$order[1], ",", spec$order[2], ")")
  } else {
    paste0("ARIMA(", spec$order[1], ",", spec$order[2], ",", spec$order[3], ")")
  }
}

resolve_innovation_spec <- function(spec_innovations, default_innovations) {
  descriptor <- spec_innovations %||% default_innovations

  if (is.function(descriptor)) {
    generator <- function(n) descriptor(n)
    theoretical <- estimate_theoretical_moments(generator)
    return(list(generator = generator,
                theoretical = theoretical,
                label = "custom",
                burn_in = .pmm2_default_burn_in))
  }

  if (!is.list(descriptor)) {
    stop("Innovation description must be a function or list.")
  }

  if (!is.null(descriptor$generator)) {
    generator <- descriptor$generator
  } else {
    generator <- innovation_generator_from_type(descriptor)
  }

  theoretical <- descriptor$theoretical
  if (is.null(theoretical)) {
    theoretical <- innovation_theoretical_moments(descriptor)
    if (any(is.na(unlist(theoretical)))) {
      theoretical <- estimate_theoretical_moments(generator)
    }
  }

  label <- descriptor$label
  if (is.null(label)) {
    label <- descriptor$type %||% "custom"
  }

  burn_in <- descriptor$burn_in %||% .pmm2_default_burn_in

  list(generator = generator,
       theoretical = theoretical,
       label = label,
       burn_in = as.integer(burn_in))
}

innovation_generator_from_type <- function(descriptor) {
  type <- match.arg(descriptor$type,
                    c("gaussian", "gamma", "student_t", "t",
                      "exponential", "chi_squared", "lognormal"))

  if (type == "gaussian") {
    return(function(n) stats::rnorm(n))
  }

  if (type %in% c("student_t", "t")) {
    df <- descriptor$df %||% 6
    if (df <= 2) {
      stop("df for student_t must be > 2")
    }
    scale <- sqrt(df / (df - 2))
    return(function(n) stats::rt(n, df = df) / scale)
  }

  if (type == "gamma") {
    shape <- descriptor$shape %||% 2
    scale <- descriptor$scale %||% 1
    return(function(n) {
      z <- stats::rgamma(n, shape = shape, scale = scale)
      centered <- z - shape * scale
      centered / (scale * sqrt(shape))
    })
  }

  if (type == "exponential") {
    rate <- descriptor$rate %||% 1
    return(function(n) {
      z <- stats::rexp(n, rate = rate)
      centered <- z - 1 / rate
      centered * rate
    })
  }

  if (type == "chi_squared") {
    df <- descriptor$df %||% 3
    return(function(n) {
      z <- stats::rchisq(n, df = df)
      centered <- z - df
      centered / sqrt(2 * df)
    })
  }

  if (type == "lognormal") {
    sigma <- descriptor$sigma %||% 0.5
    mu <- descriptor$mu %||% (-0.5 * sigma^2)
    mean_val <- exp(mu + sigma^2 / 2)
    var_val <- (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)
    sd_val <- sqrt(var_val)
    return(function(n) {
      z <- exp(stats::rnorm(n, mean = mu, sd = sigma))
      centered <- z - mean_val
      centered / sd_val
    })
  }

  stop("Unsupported innovation type: ", descriptor$type)
}

innovation_theoretical_moments <- function(descriptor) {
  type <- descriptor$type %||% "gaussian"

  if (type == "gaussian") {
    return(list(c3 = 0, c4 = 0, g = 1))
  }

  if (type %in% c("student_t", "t")) {
    df <- descriptor$df %||% 6
    if (df <= 4) {
      return(list(c3 = 0, c4 = NA_real_, g = NA_real_))
    }
    c3 <- 0
    c4 <- 6 / (df - 4)
    g <- 1 - (c3^2) / (2 + c4)
    return(list(c3 = c3, c4 = c4, g = g))
  }

  if (type == "gamma") {
    shape <- descriptor$shape %||% 2
    c3 <- 2 / sqrt(shape)
    c4 <- 6 / shape
    g <- 1 - (c3^2) / (2 + c4)
    return(list(c3 = c3, c4 = c4, g = g))
  }

  if (type == "exponential") {
    c3 <- 2
    c4 <- 6
    g <- 1 - (c3^2) / (2 + c4)
    return(list(c3 = c3, c4 = c4, g = g))
  }

  if (type == "chi_squared") {
    df <- descriptor$df %||% 3
    c3 <- sqrt(8 / df)
    c4 <- 12 / df
    g <- 1 - (c3^2) / (2 + c4)
    return(list(c3 = c3, c4 = c4, g = g))
  }

  if (type == "lognormal") {
    sigma <- descriptor$sigma %||% 0.5
    c3 <- (exp(sigma^2) + 2) * sqrt(exp(sigma^2) - 1)
    c4 <- exp(4 * sigma^2) + 2 * exp(3 * sigma^2) + 3 * exp(2 * sigma^2) - 6
    g <- 1 - (c3^2) / (2 + c4)
    return(list(c3 = c3, c4 = c4, g = g))
  }

  list(c3 = NA_real_, c4 = NA_real_, g = NA_real_)
}

estimate_theoretical_moments <- function(generator, n = 100000L) {
  sample <- generator(n)
  moments <- compute_moments(sample)
  list(c3 = moments$c3, c4 = moments$c4, g = moments$g)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
