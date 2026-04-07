## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----load_package-------------------------------------------------------------
library(EstemPMM)
set.seed(42)

## ----ar1_simulation-----------------------------------------------------------
# Model parameters
n <- 300
phi <- 0.6  # AR coefficient
intercept <- 5

# Generate skewed innovations using gamma distribution
# (Gamma(shape=2) has skewness = sqrt(2) ≈ 1.41)
innovations <- rgamma(n, shape = 2, rate = 2) - 1  # Center around 0

# Generate AR(1) series
x <- numeric(n)
x[1] <- intercept + innovations[1]
for(i in 2:n) {
  x[i] <- intercept + phi * (x[i-1] - intercept) + innovations[i]
}

# Plot the time series
plot(x, type = "l", main = "AR(1) Process with Skewed Innovations",
     ylab = "Value", xlab = "Time", col = "steelblue", lwd = 1.5)
abline(h = intercept, col = "red", lty = 2)

## ----ar1_estimation-----------------------------------------------------------
# Fit using PMM2
fit_pmm2_ar1 <- ar_pmm2(x, order = 1, method = "pmm2", include.mean = TRUE)

# Fit using CSS (classical method)
fit_css_ar1 <- ar_pmm2(x, order = 1, method = "css", include.mean = TRUE)

# Compare coefficients
cat("True AR coefficient:", phi, "\n\n")

cat("CSS Estimate:\n")
print(coef(fit_css_ar1))

cat("\nPMM2 Estimate:\n")
print(coef(fit_pmm2_ar1))

# Summary with diagnostics
summary(fit_pmm2_ar1)

## ----ar2_example--------------------------------------------------------------
# AR(2) model: x_t = 0.5*x_{t-1} - 0.3*x_{t-2} + e_t
n <- 400
phi <- c(0.5, -0.3)

# Generate t-distributed innovations (heavy tails)
innovations <- rt(n, df = 4)

# Generate AR(2) series
x <- arima.sim(n = n, model = list(ar = phi), innov = innovations)

# Fit models
fit_pmm2_ar2 <- ar_pmm2(x, order = 2, method = "pmm2", include.mean = FALSE)
fit_css_ar2 <- ar_pmm2(x, order = 2, method = "css", include.mean = FALSE)

# Compare
comparison_ar2 <- data.frame(
  Parameter = c("phi1", "phi2"),
  True = phi,
  CSS = coef(fit_css_ar2),
  PMM2 = coef(fit_pmm2_ar2)
)

print(comparison_ar2, row.names = FALSE)

## ----ma1_simulation-----------------------------------------------------------
# MA(1) model: x_t = e_t + theta*e_{t-1}
n <- 300
theta <- 0.7

# Generate chi-squared innovations (right-skewed)
innovations <- rchisq(n, df = 3) - 3

# Generate MA(1) series
x <- numeric(n)
x[1] <- innovations[1]
for(i in 2:n) {
  x[i] <- innovations[i] + theta * innovations[i-1]
}

# Fit MA(1) models
fit_pmm2_ma1 <- ma_pmm2(x, order = 1, method = "pmm2", include.mean = FALSE)
fit_css_ma1 <- ma_pmm2(x, order = 1, method = "css", include.mean = FALSE)

# Compare estimates
cat("True MA coefficient:", theta, "\n\n")

cat("CSS Estimate:\n")
print(coef(fit_css_ma1))

cat("\nPMM2 Estimate:\n")
print(coef(fit_pmm2_ma1))

# Check convergence
cat("\nPMM2 Convergence:", fit_pmm2_ma1@convergence, "\n")
cat("Iterations:", fit_pmm2_ma1@iterations, "\n")

## ----ma2_example--------------------------------------------------------------
# MA(2) model
n <- 400
theta <- c(0.6, 0.3)

# Mixed distribution: 80% normal + 20% contamination
innovations <- ifelse(runif(n) < 0.8,
                      rnorm(n),
                      rnorm(n, mean = 0, sd = 3))

# Generate MA(2) series
x <- arima.sim(n = n, model = list(ma = theta), innov = innovations)

# Fit models
fit_pmm2_ma2 <- ma_pmm2(x, order = 2, method = "pmm2", include.mean = FALSE)
fit_css_ma2 <- ma_pmm2(x, order = 2, method = "css", include.mean = FALSE)

# Comparison
comparison_ma2 <- data.frame(
  Parameter = c("theta1", "theta2"),
  True = theta,
  CSS = coef(fit_css_ma2),
  PMM2 = coef(fit_pmm2_ma2)
)

print(comparison_ma2, row.names = FALSE)

## ----arma11_simulation--------------------------------------------------------
# ARMA(1,1) model
n <- 500
phi <- 0.5    # AR coefficient
theta <- 0.4  # MA coefficient

# Generate exponentially distributed innovations (right-skewed)
innovations <- rexp(n, rate = 1) - 1

# Generate ARMA(1,1) series
x <- arima.sim(n = n,
               model = list(ar = phi, ma = theta),
               innov = innovations)

# Fit models
fit_pmm2_arma <- arma_pmm2(x, order = c(1, 1), method = "pmm2", include.mean = FALSE)
fit_css_arma <- arma_pmm2(x, order = c(1, 1), method = "css", include.mean = FALSE)

# Extract coefficients
cat("True parameters: AR =", phi, ", MA =", theta, "\n\n")

cat("CSS Estimates:\n")
print(coef(fit_css_arma))

cat("\nPMM2 Estimates:\n")
print(coef(fit_pmm2_arma))

# Summary
summary(fit_pmm2_arma)

## ----arma_diagnostics, fig.width=7, fig.height=6------------------------------
# Plot residuals and ACF
plot(fit_pmm2_arma)
title("ARMA(1,1) Model Diagnostics")

## ----arima111_simulation------------------------------------------------------
# Generate non-stationary series requiring differencing
n <- 400
phi <- 0.4
theta <- 0.3

# Skewed innovations
innovations <- rgamma(n, shape = 3, rate = 3) - 1

# Generate integrated series: ARIMA(1,1,1)
# First generate stationary ARMA(1,1)
x_stationary <- arima.sim(n = n,
                          model = list(ar = phi, ma = theta),
                          innov = innovations)

# Integrate once (cumulative sum)
x <- cumsum(x_stationary)

# Plot the non-stationary series
plot(x, type = "l", main = "ARIMA(1,1,1) Process",
     ylab = "Value", xlab = "Time", col = "darkgreen", lwd = 1.5)

## ----arima_estimation---------------------------------------------------------
# Fit ARIMA(1,1,1) using PMM2
fit_pmm2_arima <- arima_pmm2(x, order = c(1, 1, 1), method = "pmm2", include.mean = FALSE)

# Fit using classical CSS method
fit_css_arima <- arima_pmm2(x, order = c(1, 1, 1), method = "css", include.mean = FALSE)

# Compare coefficients
cat("True parameters: AR =", phi, ", MA =", theta, "\n\n")

cat("CSS Estimates:\n")
print(coef(fit_css_arima))

cat("\nPMM2 Estimates:\n")
print(coef(fit_pmm2_arima))

# Summary
summary(fit_pmm2_arima)

## ----ts_bootstrap_inference---------------------------------------------------
# Perform bootstrap inference for the ARMA(1,1) model
boot_results <- ts_pmm2_inference(fit_pmm2_arma,
                                  B = 500,
                                  seed = 123,
                                  method = "block",
                                  block_length = 20,
                                  parallel = FALSE)

# Display summary with confidence intervals
summary(boot_results)

## ----monte_carlo_comparison, eval=FALSE---------------------------------------
# # Define model specifications
# specs <- list(
#   list(
#     model = "ma",
#     order = 1,
#     theta = 0.7,
#     label = "MA(1)",
#     innovations = list(type = "gamma", shape = 2)
#   ),
#   list(
#     model = "arma",
#     order = c(1, 1),
#     theta = list(ar = 0.5, ma = 0.4),
#     label = "ARMA(1,1)",
#     innovations = list(type = "student_t", df = 4)
#   )
# )
# 
# # Run Monte Carlo comparison
# mc_results <- pmm2_monte_carlo_compare(
#   model_specs = specs,
#   methods = c("css", "pmm2"),
#   n = 300,
#   n_sim = 1000,
#   seed = 456,
#   progress = FALSE
# )
# 
# # Display results
# print(mc_results$summary)
# print(mc_results$gain)

## ----prediction---------------------------------------------------------------
# Forecast 10 steps ahead for the ARMA(1,1) model
forecast_horizon <- 10

pred_raw <- predict(fit_pmm2_arma, n.ahead = forecast_horizon)
predictions <- if (is.list(pred_raw)) as.numeric(pred_raw$pred) else as.numeric(pred_raw)

# Display forecasts
cat("Forecasts for next", forecast_horizon, "periods:\n")
print(predictions)

# Plot original series with forecasts
n_plot <- 100
plot_data <- tail(x, n_plot)
plot(seq_along(plot_data), plot_data, type = "l",
     xlim = c(1, n_plot + forecast_horizon),
     ylim = range(c(plot_data, predictions)),
     main = "ARMA(1,1) Forecasts",
     xlab = "Time", ylab = "Value",
     col = "steelblue", lwd = 1.5)

# Add forecasts
lines(n_plot + seq_len(forecast_horizon), predictions,
      col = "red", lwd = 2, lty = 2)

legend("topleft",
       legend = c("Observed", "Forecast"),
       col = c("steelblue", "red"),
       lty = c(1, 2), lwd = c(1.5, 2))

## ----custom_innovations, eval=FALSE-------------------------------------------
# # Custom innovation generator
# my_innovations <- function(n) {
#   # Mixture: 90% Gaussian + 10% extreme outliers
#   base <- rnorm(n)
#   outliers <- sample(c(-5, 5), n, replace = TRUE, prob = c(0.05, 0.05))
#   ifelse(runif(n) < 0.9, base, outliers)
# }
# 
# # Use in Monte Carlo
# specs <- list(
#   list(
#     model = "ar",
#     order = 1,
#     theta = 0.6,
#     innovations = list(generator = my_innovations)
#   )
# )
# 
# mc_results <- pmm2_monte_carlo_compare(specs, methods = c("css", "pmm2"),
#                                        n = 300, n_sim = 500)

## ----sar_example--------------------------------------------------------------
# Simulate SAR(1,1)_12 with seasonal period = 12
n <- 300
phi <- 0.6       # Non-seasonal AR coefficient
Phi <- 0.4       # Seasonal AR coefficient (lag 12)
s <- 12          # Seasonal period

# Generate gamma-distributed innovations (right-skewed)
innovations <- rgamma(n, shape = 2, rate = 2) - 1

# Generate SAR series
y <- numeric(n)
for (t in (s+2):n) {
  y[t] <- phi * y[t-1] + Phi * y[t-s] + innovations[t]
}

# Plot the series
plot(y, type = "l", main = "SAR(1,1)_12 Process with Skewed Innovations",
     ylab = "Value", xlab = "Time", col = "steelblue", lwd = 1.5)

# Fit SAR model with PMM2
fit_sar_pmm2 <- sar_pmm2(
  y,
  order = c(1, 1),
  season = list(period = 12),
  method = "pmm2"
)

# Compare with OLS
fit_sar_ols <- sar_pmm2(
  y,
  order = c(1, 1),
  season = list(period = 12),
  method = "ols"
)

# Display results
cat("True parameters: phi =", phi, ", Phi =", Phi, "\n\n")

cat("OLS Estimates:\n")
print(coef(fit_sar_ols))

cat("\nPMM2 Estimates:\n")
print(coef(fit_sar_pmm2))

# Summary with diagnostics
summary(fit_sar_pmm2)

## ----sma_example--------------------------------------------------------------
# Simulate SMA(1)_4 (quarterly seasonal pattern)
n <- 200
Theta <- 0.6     # Seasonal MA coefficient
s <- 4           # Quarterly data

# Generate exponentially distributed innovations
innovations <- rexp(n, rate = 1) - 1

# Generate SMA series
y <- numeric(n)
for (t in 1:s) {
  y[t] <- innovations[t]
}
for (t in (s+1):n) {
  y[t] <- innovations[t] + Theta * innovations[t-s]
}

# Plot the series
plot(y, type = "l", main = "SMA(1)_4 Process",
     ylab = "Value", xlab = "Time", col = "darkgreen", lwd = 1.5)

# Fit SMA model with PMM2
fit_sma_pmm2 <- sma_pmm2(
  y,
  order = 1,
  season = list(period = 4),
  method = "pmm2"
)

# Fit with CSS for comparison
fit_sma_css <- sma_pmm2(
  y,
  order = 1,
  season = list(period = 4),
  method = "css"
)

# Display results
cat("True parameter: Theta =", Theta, "\n\n")

cat("CSS Estimate:\n")
print(coef(fit_sma_css))

cat("\nPMM2 Estimate:\n")
print(coef(fit_sma_pmm2))

cat("\nConvergence:", fit_sma_pmm2@convergence, "\n")
cat("Iterations:", fit_sma_pmm2@iterations, "\n")

## ----compare_seasonal, eval=FALSE---------------------------------------------
# # Systematic comparison of SAR estimation methods
# compare_sar_methods(y, order = c(1, 1), period = 12)
# 
# # Universal comparison for non-seasonal models
# compare_ts_methods(
#   y,
#   model_type = "arima",
#   order = c(1, 0, 1)
# )

