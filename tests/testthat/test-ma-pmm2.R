test_that("estpmm_style_ma works for MA(1) models", {
    set.seed(123)
    n <- 200
    theta <- 0.6

    # Generate MA(1) with exponential errors (asymmetric)
    innovations <- rexp(n, rate = 1) - 1
    x <- arima.sim(n = n, list(ma = theta), innov = innovations)

    # Fit using new PMM2 method
    fit <- EstemPMM:::estpmm_style_ma(x, q = 1, include.mean = FALSE, verbose = FALSE)

    expect_type(fit, "list")
    expect_true(fit$convergence)
    expect_length(fit$ma_coef, 1)
    expect_equal(length(fit$innovations), n)

    # Check that estimates are reasonable
    expect_equal(fit$ma_coef[[1]], theta, tolerance = 0.2)

    # Check method string
    expect_equal(fit$method, "EstemPMM-style PMM2")
})

test_that("estpmm_style_sma works for SMA(1) models", {
    set.seed(456)
    n <- 200
    Theta <- 0.6
    s <- 4

    # Generate SMA(1)_4
    innovations <- rexp(n, rate = 1) - 1
    # Manual generation for SMA
    x <- numeric(n)
    for (t in 1:n) {
        sma_term <- if (t > s) Theta * innovations[t - s] else 0
        x[t] <- innovations[t] + sma_term
    }

    # Fit using new PMM2 method
    fit <- EstemPMM:::estpmm_style_sma(x, Q = 1, s = s, include.mean = FALSE, verbose = FALSE)

    expect_type(fit, "list")
    expect_true(fit$convergence)
    expect_length(fit$sma_coef, 1)

    # Check estimates
    expect_equal(fit$sma_coef[[1]], Theta, tolerance = 0.2)
})

test_that("sarima_pmm2 integrates ma_method='pmm2' correctly", {
    set.seed(789)
    n <- 200

    # MA(1) case
    innovations <- rexp(n, rate = 1) - 1
    x <- arima.sim(n = n, list(ma = 0.5), innov = innovations)

    # Fit with ma_method = "pmm2"
    fit <- sarima_pmm2(x,
        order = c(0, 0, 1, 0), seasonal = list(order = c(0, 0), period = 1),
        ma_method = "pmm2", include.mean = FALSE
    )

    expect_s4_class(fit, "SARIMAPMM2")
    expect_true(fit@convergence)
    expect_equal(length(fit@coefficients), 1)

    # Fit with ma_method = "mle" (default)
    fit_mle <- sarima_pmm2(x,
        order = c(0, 0, 1, 0), seasonal = list(order = c(0, 0), period = 1),
        ma_method = "mle", include.mean = FALSE
    )

    expect_s4_class(fit_mle, "SARIMAPMM2")

    # Ensure results are slightly different (since methods differ)
    # Note: For normal data they might be close, but for exp data they should differ
    expect_false(isTRUE(all.equal(fit@coefficients, fit_mle@coefficients)))
})

test_that("estpmm_style_ma_sma works for MA(1)+SMA(1) models", {
    set.seed(999)
    n <- 200
    theta_ma <- 0.5
    theta_sma <- 0.6
    s <- 12

    # Generate MA(1)+SMA(1)_12 with asymmetric errors
    innovations <- rexp(n, rate = 1) - 1
    x <- numeric(n)

    for (t in 1:n) {
        ma_term <- if (t > 1) theta_ma * innovations[t - 1] else 0
        sma_term <- if (t > s) theta_sma * innovations[t - s] else 0
        x[t] <- innovations[t] + ma_term + sma_term
    }

    # Fit using new MA+SMA PMM2 method
    fit <- EstemPMM:::estpmm_style_ma_sma(x, q = 1, Q = 1, s = s, include.mean = FALSE, verbose = FALSE)

    expect_type(fit, "list")
    expect_true(fit$convergence)
    expect_length(fit$ma_coef, 1)
    expect_length(fit$sma_coef, 1)
    expect_equal(length(fit$innovations), n)

    # Check that estimates are reasonable (relaxed tolerance given complexity)
    expect_equal(fit$ma_coef[[1]], theta_ma, tolerance = 0.4)
    expect_equal(fit$sma_coef[[1]], theta_sma, tolerance = 0.4)

    # Check method string
    expect_equal(fit$method, "EstemPMM-style PMM2")
})

test_that("sarima_pmm2 integrates MA+SMA with ma_method='pmm2'", {
    set.seed(1111)
    n <- 200
    theta_ma <- 0.4
    theta_sma <- 0.5
    s <- 4

    # Generate MA(1)+SMA(1)_4
    innovations <- rgamma(n, shape = 2, scale = 1) - 2 # asymmetric
    x <- numeric(n)

    for (t in 1:n) {
        ma_term <- if (t > 1) theta_ma * innovations[t - 1] else 0
        sma_term <- if (t > s) theta_sma * innovations[t - s] else 0
        x[t] <- innovations[t] + ma_term + sma_term
    }

    # Fit with ma_method = "pmm2"
    fit_pmm2 <- sarima_pmm2(x,
        order = c(0, 0, 1, 1), seasonal = list(order = c(0, 0), period = s),
        ma_method = "pmm2", include.mean = FALSE
    )

    expect_s4_class(fit_pmm2, "SARIMAPMM2")
    expect_true(fit_pmm2@convergence)
    expect_equal(length(fit_pmm2@coefficients), 2) # ma1 + sma1
    expect_named(fit_pmm2@coefficients, c("ma1", "sma1"))

    # Fit with ma_method = "mle" (default)
    fit_mle <- sarima_pmm2(x,
        order = c(0, 0, 1, 1), seasonal = list(order = c(0, 0), period = s),
        ma_method = "mle", include.mean = FALSE
    )

    expect_s4_class(fit_mle, "SARIMAPMM2")

    # Both methods should converge successfully
    # Note: Differences may be small with gamma distribution; main check is convergence
    expect_true(fit_mle@convergence)
    expect_equal(length(fit_mle@coefficients), 2)
})

test_that("MA+SMA PMM2 handles different orders correctly", {
    set.seed(2222)
    n <- 150
    s <- 12

    # MA(2)+SMA(1)
    innovations <- rexp(n, rate = 1) - 1
    x <- numeric(n)
    theta_ma <- c(0.3, 0.2)
    theta_sma <- c(0.5)

    for (t in 1:n) {
        ma_term <- 0
        if (t > 1) ma_term <- ma_term + theta_ma[1] * innovations[t - 1]
        if (t > 2) ma_term <- ma_term + theta_ma[2] * innovations[t - 2]
        sma_term <- if (t > s) theta_sma[1] * innovations[t - s] else 0
        x[t] <- innovations[t] + ma_term + sma_term
    }

    # Fit MA(2)+SMA(1)
    fit <- EstemPMM:::estpmm_style_ma_sma(x, q = 2, Q = 1, s = s, include.mean = FALSE, verbose = FALSE)

    expect_length(fit$ma_coef, 2)
    expect_length(fit$sma_coef, 1)
    expect_true(fit$convergence)

    # Coefficients should be in ballpark
    expect_equal(fit$ma_coef[1], theta_ma[1], tolerance = 0.3)
    expect_equal(fit$sma_coef[1], theta_sma[1], tolerance = 0.3)
})
