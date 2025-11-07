set.seed(42)
n <- 80
x <- rnorm(n)
y <- 2 + 1.5 * x + rgamma(n, shape = 2, scale = 1) - 2
test_data <- data.frame(x = x, y = y)

test_that("pmm2_inference returns bootstrap summary", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  inf <- pmm2_inference(fit, formula = y ~ x, data = test_data, B = 20, seed = 123)

  expect_s3_class(inf, "data.frame")
  expect_true(all(c("Estimate", "Std.Error", "t.value", "p.value") %in% names(inf)))
  expect_true(nrow(inf) >= 2)
  expect_true(all(is.finite(inf$Std.Error)))
})

test_that("pmm2_inference confidence intervals are ordered", {
  fit <- lm_pmm2(y ~ x, data = test_data)

  inf <- pmm2_inference(fit, formula = y ~ x, data = test_data, B = 18, seed = 321)

  expect_true(all(inf$conf.low < inf$conf.high))
})

test_that("ts_pmm2_inference works for AR models", {
  ts_data <- as.numeric(arima.sim(model = list(ar = 0.6), n = 120))
  fit <- ar_pmm2(ts_data, order = 1)

  inf <- ts_pmm2_inference(fit, B = 15, seed = 99)

  expect_s3_class(inf, "data.frame")
  expect_true(all(c("Estimate", "Std.Error", "conf.low", "conf.high") %in% names(inf)))
  expect_gt(nrow(inf), 0)
})
