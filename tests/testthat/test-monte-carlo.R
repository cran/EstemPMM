test_that("pmm2_monte_carlo_compare returns structured results", {
  specs <- list(
    list(model = "ma", order = 1, theta = 0.6, label = "MA(1)")
  )

  res <- pmm2_monte_carlo_compare(
    model_specs = specs,
    methods = c("css", "pmm2"),
    n = 40,
    n_sim = 4,
    seed = 321,
    progress = FALSE
  )

  expect_type(res, "list")
  expect_true(all(c("parameter_results", "summary", "gain") %in% names(res)))

  param_df <- res$parameter_results
  expect_s3_class(param_df, "data.frame")
  expect_true(all(c("model", "method", "parameter", "mse", "mse_ratio",
                    "theoretical_c3", "empirical_c3") %in% names(param_df)))
  expect_true(any(param_df$method == "css"))
  expect_true(any(param_df$method == "pmm2"))

  gain_df <- res$gain
  expect_s3_class(gain_df, "data.frame")
  expect_true(all(c("model", "theoretical_g", "empirical_g", "observed_ratio") %in% names(gain_df)))
})

test_that("gamma innovations produce expected theoretical moments", {
  shape_val <- 4
  specs <- list(
    list(
      model = "ma",
      order = 1,
      theta = 0.5,
      label = "MA(1)",
      innovations = list(type = "gamma", shape = shape_val)
    )
  )

  res <- pmm2_monte_carlo_compare(
    model_specs = specs,
    methods = c("css", "pmm2"),
    n = 30,
    n_sim = 3,
    seed = 42,
    progress = FALSE
  )

  gain_df <- res$gain
  expect_true(nrow(gain_df) == 1)

  theo_c3 <- gain_df$theoretical_c3
  theo_c4 <- gain_df$theoretical_c4
  expect_equal(theo_c3, 2 / sqrt(shape_val), tolerance = 1e-8)
  expect_equal(theo_c4, 6 / shape_val, tolerance = 1e-8)
  expect_equal(gain_df$theoretical_g, 1 - (theo_c3^2) / (2 + theo_c4), tolerance = 1e-8)
})
