## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----load_library-------------------------------------------------------------
library(EstemPMM)
library(ggplot2)
set.seed(2025)

## ----sar_1_1_12---------------------------------------------------------------
# Simulate SAR(1,1)_12 model
n <- 360  # 30 years of monthly data
phi <- 0.5      # Non-seasonal AR(1)
Phi <- 0.6      # Seasonal AR(1) at lag 12
s <- 12         # Monthly seasonal period

# Generate gamma-distributed innovations (right-skewed)
# This mimics asymmetric economic shocks
innovations <- rgamma(n + 100, shape = 2, scale = 1) - 2

# Generate SAR series
y <- numeric(n + 100)
for (t in (s + 2):(n + 100)) {
  y[t] <- phi * y[t-1] + Phi * y[t-s] + innovations[t]
}
y <- y[101:(n + 100)]  # Remove burn-in

# Convert to time series object for better plotting
y_ts <- ts(y, frequency = 12, start = c(2000, 1))

# Plot the series
plot(y_ts, main = "SAR(1,1)_12 Process with Skewed Innovations",
     ylab = "Value", xlab = "Year", col = "steelblue", lwd = 1.5)
abline(h = 0, col = "gray", lty = 2)

## ----sar_seasonal_pattern-----------------------------------------------------
# Extract seasonal patterns
months <- rep(1:12, length.out = length(y))
seasonal_df <- data.frame(
  Month = factor(months, labels = month.abb),
  Value = y
)

# Boxplot by month
boxplot(Value ~ Month, data = seasonal_df,
        main = "Seasonal Pattern in SAR(1,1)_12 Data",
        xlab = "Month", ylab = "Value",
        col = "lightblue", border = "steelblue")

## ----fit_sar------------------------------------------------------------------
# Fit SAR(1,1)_12 with PMM2
fit_sar_pmm2 <- sar_pmm2(
  y,
  order = c(1, 1),           # (non-seasonal AR order, seasonal AR order)
  season = list(period = 12), # Monthly seasonality
  method = "pmm2",
  verbose = FALSE
)

# Fit with OLS for comparison
fit_sar_ols <- sar_pmm2(
  y,
  order = c(1, 1),
  season = list(period = 12),
  method = "ols",
  verbose = FALSE
)

# Fit with CSS for comparison
fit_sar_css <- sar_pmm2(
  y,
  order = c(1, 1),
  season = list(period = 12),
  method = "css",
  verbose = FALSE
)

# Compare estimates
cat("True parameters:\n")
cat("  phi (AR1):", phi, "\n")
cat("  Phi (SAR1):", Phi, "\n\n")

cat("OLS Estimates:\n")
print(coef(fit_sar_ols))

cat("\nCSS Estimates:\n")
print(coef(fit_sar_css))

cat("\nPMM2 Estimates:\n")
print(coef(fit_sar_pmm2))

# Display PMM2 summary with diagnostics
cat("\n")
summary(fit_sar_pmm2)

## ----sar_g_factor-------------------------------------------------------------
# Extract moments from fit
m2 <- fit_sar_pmm2@m2
m3 <- fit_sar_pmm2@m3
m4 <- fit_sar_pmm2@m4

# Compute g-factor
g_result <- pmm2_variance_factor(m2, m3, m4)
cat("PMM2 variance correction factor:\n")
cat("  g =", round(g_result$g, 4), "\n")
cat("  Theoretical variance reduction:", 
    round((1 - g_result$g) * 100, 2), "%\n")

# For Gaussian innovations, g ≈ 1 (no improvement)
# For skewed/heavy-tailed innovations, g < 1 (PMM2 more efficient)

## ----sar_diagnostics, fig.height=7--------------------------------------------
# Create diagnostic plots
par(mfrow = c(2, 2))

# 1. Residuals over time
residuals_pmm2 <- fit_sar_pmm2@residuals
plot(residuals_pmm2, type = "l",
     main = "PMM2 Residuals",
     xlab = "Time", ylab = "Residual",
     col = "steelblue")
abline(h = 0, col = "red", lty = 2)

# 2. ACF of residuals
acf(residuals_pmm2, main = "ACF of PMM2 Residuals", 
    col = "steelblue", lwd = 2)

# 3. Q-Q plot
qqnorm(residuals_pmm2, main = "Q-Q Plot of PMM2 Residuals",
       pch = 20, col = "steelblue")
qqline(residuals_pmm2, col = "red", lwd = 2)

# 4. Histogram of residuals
hist(residuals_pmm2, breaks = 30, col = "lightblue",
     border = "steelblue", main = "Distribution of PMM2 Residuals",
     xlab = "Residual", freq = FALSE)
lines(density(residuals_pmm2), col = "darkblue", lwd = 2)

par(mfrow = c(1, 1))

## ----sma_1_4------------------------------------------------------------------
# Simulate SMA(1)_4 model (quarterly seasonality)
n <- 200
Theta <- 0.6    # Seasonal MA(1) coefficient
s <- 4          # Quarterly period

# Generate exponentially distributed innovations (right-skewed)
innovations <- rexp(n + 50, rate = 1) - 1

# Generate SMA series
y <- numeric(n + 50)
for (t in 1:s) {
  y[t] <- innovations[t]
}
for (t in (s + 1):(n + 50)) {
  y[t] <- innovations[t] + Theta * innovations[t - s]
}
y <- y[51:(n + 50)]  # Remove burn-in

# Convert to quarterly time series
y_ts <- ts(y, frequency = 4, start = c(2000, 1))

# Plot
plot(y_ts, main = "SMA(1)_4 Process with Exponential Innovations",
     ylab = "Value", xlab = "Year", col = "darkgreen", lwd = 1.5)
abline(h = 0, col = "gray", lty = 2)

## ----sma_quarterly_pattern----------------------------------------------------
quarters <- rep(1:4, length.out = length(y))
seasonal_df_q <- data.frame(
  Quarter = factor(quarters, labels = paste("Q", 1:4, sep = "")),
  Value = y
)

boxplot(Value ~ Quarter, data = seasonal_df_q,
        main = "Seasonal Pattern in SMA(1)_4 Data",
        xlab = "Quarter", ylab = "Value",
        col = "lightgreen", border = "darkgreen")

## ----fit_sma------------------------------------------------------------------
# Fit SMA(1)_4 with PMM2
fit_sma_pmm2 <- sma_pmm2(
  y,
  order = 1,                 # Seasonal MA order
  season = list(period = 4), # Quarterly seasonality
  method = "pmm2",
  verbose = FALSE
)

# Fit with CSS for comparison
fit_sma_css <- sma_pmm2(
  y,
  order = 1,
  season = list(period = 4),
  method = "css",
  verbose = FALSE
)

# Compare estimates
cat("True parameter: Theta =", Theta, "\n\n")

cat("CSS Estimate:\n")
print(coef(fit_sma_css))

cat("\nPMM2 Estimate:\n")
print(coef(fit_sma_pmm2))

cat("\n")
summary(fit_sma_pmm2)

## ----sma_convergence----------------------------------------------------------
cat("PMM2 Convergence:\n")
cat("  Converged:", fit_sma_pmm2@convergence, "\n")
cat("  Iterations:", fit_sma_pmm2@iterations, "\n")
cat("  Residual variance (m2):", round(fit_sma_pmm2@m2, 4), "\n")

# Extract g-factor
g_sma <- pmm2_variance_factor(fit_sma_pmm2@m2, 
                              fit_sma_pmm2@m3, 
                              fit_sma_pmm2@m4)
cat("\nEfficiency:\n")
cat("  g-factor:", round(g_sma$g, 4), "\n")
cat("  Variance reduction:", round((1 - g_sma$g) * 100, 2), "%\n")

## ----sar_monte_carlo, eval=FALSE----------------------------------------------
# # This demonstrates the Monte Carlo approach
# # (Set eval=TRUE to run, but it takes several minutes)
# 
# n_sims <- 500
# n <- 200
# 
# sar_results <- data.frame(
#   method = character(),
#   phi_estimate = numeric(),
#   Phi_estimate = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# for (i in 1:n_sims) {
#   # Generate data
#   innov <- rgamma(n + 100, shape = 2, scale = 1) - 2
#   y_sim <- numeric(n + 100)
#   for (t in 13:(n + 100)) {
#     y_sim[t] <- 0.5 * y_sim[t-1] + 0.6 * y_sim[t-12] + innov[t]
#   }
#   y_sim <- y_sim[101:(n + 100)]
# 
#   # Fit PMM2
#   fit_pmm2 <- sar_pmm2(y_sim, order = c(1, 1),
#                        season = list(period = 12),
#                        method = "pmm2", verbose = FALSE)
# 
#   # Fit CSS
#   fit_css <- sar_pmm2(y_sim, order = c(1, 1),
#                       season = list(period = 12),
#                       method = "css", verbose = FALSE)
# 
#   # Store results
#   sar_results <- rbind(sar_results,
#     data.frame(method = "PMM2",
#                phi_estimate = coef(fit_pmm2)[1],
#                Phi_estimate = coef(fit_pmm2)[2]),
#     data.frame(method = "CSS",
#                phi_estimate = coef(fit_css)[1],
#                Phi_estimate = coef(fit_css)[2])
#   )
# }
# 
# # Compute variance reduction
# var_pmm2_phi <- var(sar_results$phi_estimate[sar_results$method == "PMM2"])
# var_css_phi <- var(sar_results$phi_estimate[sar_results$method == "CSS"])
# var_reduction_phi <- (1 - var_pmm2_phi / var_css_phi) * 100
# 
# var_pmm2_Phi <- var(sar_results$Phi_estimate[sar_results$method == "PMM2"])
# var_css_Phi <- var(sar_results$Phi_estimate[sar_results$method == "CSS"])
# var_reduction_Phi <- (1 - var_pmm2_Phi / var_css_Phi) * 100
# 
# cat("SAR(1,1)_12 Monte Carlo Results (n=200, 500 sims):\n")
# cat("  Variance reduction for phi:", round(var_reduction_phi, 2), "%\n")
# cat("  Variance reduction for Phi:", round(var_reduction_Phi, 2), "%\n")

## ----sarma_example------------------------------------------------------------
# Simulate SARMA(1,0,0,1)_12 model
# This is SAR(1) × SMA(1) with period 12
n <- 300
phi <- 0.5
Theta <- 0.4
s <- 12

# Generate lognormal innovations (right-skewed)
innovations <- rlnorm(n + 100, meanlog = 0, sdlog = 0.5) - exp(0.125)

# For simplicity, use arima.sim
y <- arima.sim(
  n = n,
  model = list(
    ar = phi,
    seasonal = list(order = c(1, 0, 1), period = 12, 
                    sar = 0, sma = Theta)
  ),
  innov = innovations[1:n]
)

# Fit SARMA model with PMM2
fit_sarma_pmm2 <- sarma_pmm2(
  y,
  order = c(1, 1, 0, 1),     # (p, P, q, Q)
  season = list(period = 12),
  method = "pmm2",
  verbose = FALSE
)

# Fit with CSS
fit_sarma_css <- sarma_pmm2(
  y,
  order = c(1, 1, 0, 1),
  season = list(period = 12),
  method = "css",
  verbose = FALSE
)

# Compare
cat("True parameters: phi =", phi, ", Theta =", Theta, "\n\n")

cat("CSS Estimates:\n")
print(coef(fit_sarma_css))

cat("\nPMM2 Estimates:\n")
print(coef(fit_sarma_pmm2))

## ----sarima_example, eval=FALSE-----------------------------------------------
# # Simulate SARIMA(1,1,0)×(1,1,1)_12
# # This includes regular differencing (d=1) and seasonal differencing (D=1)
# n <- 400
# phi <- 0.4
# Phi <- 0.5
# Theta <- 0.6
# 
# # Generate series (using standard arima.sim)
# innovations <- rgamma(n + 100, shape = 2, scale = 1) - 2
# 
# y <- arima.sim(
#   n = n,
#   model = list(
#     order = c(1, 1, 0),
#     ar = phi,
#     seasonal = list(order = c(1, 1, 1), period = 12,
#                     sar = Phi, sma = Theta)
#   ),
#   innov = innovations[1:n]
# )
# 
# # Plot the non-stationary series
# plot(y, type = "l", main = "SARIMA(1,1,0)×(1,1,1)_12 Process",
#      ylab = "Value", xlab = "Time", col = "purple", lwd = 1.5)
# 
# # Fit SARIMA with PMM2
# fit_sarima_pmm2 <- sarima_pmm2(
#   y,
#   order = c(1, 1, 0, 1),       # Model orders: (p, P, q, Q)
#   seasonal = list(
#     order = c(1, 1),           # Differencing orders: (d, D)
#     period = 12
#   ),
#   method = "pmm2",
#   verbose = FALSE
# )
# 
# # Fit with CSS
# fit_sarima_css <- sarima_pmm2(
#   y,
#   order = c(1, 1, 0, 1),       # Model orders: (p, P, q, Q)
#   seasonal = list(order = c(1, 1), period = 12),
#   method = "css",
#   verbose = FALSE
# )
# 
# # Compare
# cat("True parameters:\n")
# cat("  phi =", phi, ", Phi =", Phi, ", Theta =", Theta, "\n\n")
# 
# cat("CSS Estimates:\n")
# print(coef(fit_sarima_css))
# 
# cat("\nPMM2 Estimates:\n")
# print(coef(fit_sarima_pmm2))
# 
# summary(fit_sarima_pmm2)

## ----method_comparison--------------------------------------------------------
# For the SAR(1,1)_12 data, compare all available methods
methods <- c("ols", "css", "pmm2")

comparison_results <- data.frame(
  Method = methods,
  phi_estimate = numeric(length(methods)),
  Phi_estimate = numeric(length(methods)),
  AIC = numeric(length(methods))
)

for (i in seq_along(methods)) {
  fit <- sar_pmm2(y, order = c(1, 1), season = list(period = 12),
                 method = methods[i], verbose = FALSE)
  comparison_results$phi_estimate[i] <- coef(fit)[1]
  comparison_results$Phi_estimate[i] <- coef(fit)[2]
  comparison_results$AIC[i] <- AIC(fit)
}

print(comparison_results)

# Best method by AIC
best_method <- comparison_results$Method[which.min(comparison_results$AIC)]
cat("\nBest method by AIC:", best_method, "\n")

## ----information_criteria-----------------------------------------------------
# Compare different SAR orders
models <- list(
  list(name = "SAR(1,0)", order = c(1, 0)),
  list(name = "SAR(0,1)", order = c(0, 1)),
  list(name = "SAR(1,1)", order = c(1, 1)),
  list(name = "SAR(2,1)", order = c(2, 1))
)

ic_results <- data.frame(
  Model = character(),
  AIC = numeric(),
  BIC = numeric(),
  stringsAsFactors = FALSE
)

for (model in models) {
  fit <- sar_pmm2(y, order = model$order, 
                 season = list(period = 12),
                 method = "pmm2", verbose = FALSE)
  
  ic_results <- rbind(ic_results, data.frame(
    Model = model$name,
    AIC = AIC(fit),
    BIC = BIC(fit)
  ))
}

print(ic_results)

# Best model
best_aic <- ic_results$Model[which.min(ic_results$AIC)]
best_bic <- ic_results$Model[which.min(ic_results$BIC)]

cat("\nBest model by AIC:", best_aic, "\n")
cat("Best model by BIC:", best_bic, "\n")

## ----airline_data-------------------------------------------------------------
# Use classic airline passengers dataset
# (In practice, load with: data(AirPassengers))
# Here we'll simulate similar data

# Simulate monthly airline passengers with trend and seasonality
n <- 144  # 12 years of monthly data
t <- 1:n

# Trend component
trend <- 100 + 2 * t

# Seasonal component (annual cycle)
seasonal <- 20 * sin(2 * pi * (t - 1) / 12)

# Random component (log-normal innovations for multiplicative effect)
set.seed(123)
random <- rlnorm(n, meanlog = 0, sdlog = 0.1)

# Combine (multiplicative model)
passengers <- (trend + seasonal) * random

# Convert to time series
passengers_ts <- ts(passengers, frequency = 12, start = c(2000, 1))

# Plot
plot(passengers_ts, main = "Simulated Airline Passengers",
     ylab = "Passengers", xlab = "Year",
     col = "darkblue", lwd = 1.5)

# Decompose
decomp <- decompose(passengers_ts)
plot(decomp)

## ----airline_modeling---------------------------------------------------------
# Take log to stabilize variance
log_passengers <- log(passengers_ts)

# First differencing for trend
diff_passengers <- diff(log_passengers)

# Seasonal differencing
diff_seasonal <- diff(diff_passengers, lag = 12)

# Plot differenced series
plot(diff_seasonal, main = "Differenced Log Passengers",
     ylab = "Differenced Value", xlab = "Year",
     col = "steelblue", lwd = 1.5)
abline(h = 0, col = "red", lty = 2)

# Check ACF and PACF
par(mfrow = c(1, 2))
acf(diff_seasonal, main = "ACF of Differenced Series", na.action = na.pass)
pacf(diff_seasonal, main = "PACF of Differenced Series", na.action = na.pass)
par(mfrow = c(1, 1))

# Fit SARIMA(0,1,1)×(0,1,1)_12 (classic airline model)
fit_airline_pmm2 <- sarima_pmm2(
  log_passengers,
  order = c(0, 0, 1, 1),           # Model orders: (p, P, q, Q)
  seasonal = list(
    order = c(1, 1),               # Differencing orders: (d, D)
    period = 12
  ),
  method = "pmm2",
  verbose = FALSE
)

# Summary
summary(fit_airline_pmm2)

# Forecast
forecast_horizon <- 24  # 2 years ahead
forecast_result <- predict(fit_airline_pmm2, n.ahead = forecast_horizon)

# Extract predictions
if (is.list(forecast_result)) {
  forecasts <- as.numeric(forecast_result$pred)
} else {
  forecasts <- as.numeric(forecast_result)
}

# Plot with forecasts
plot(log_passengers, xlim = c(2000, 2014),
     ylim = range(c(log_passengers, forecasts)),
     main = "Log Passengers with 2-Year Forecast",
     ylab = "Log(Passengers)", xlab = "Year",
     col = "darkblue", lwd = 1.5)
lines(ts(forecasts, frequency = 12, start = c(2012, 1)),
      col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Observed", "Forecast"),
       col = c("darkblue", "red"), lty = c(1, 2), lwd = c(1.5, 2))

## ----diagnostics_checklist----------------------------------------------------
# Example diagnostic workflow
residuals_check <- fit_sar_pmm2@residuals

# 1. Visual inspection
plot(residuals_check, type = "l", main = "Residuals Over Time",
     col = "steelblue", ylab = "Residual", xlab = "Time")
abline(h = 0, col = "red", lty = 2)

# 2. ACF check
acf(residuals_check, main = "ACF of Residuals", lag.max = 48)

# 3. Ljung-Box test
lb_test <- Box.test(residuals_check, lag = 20, type = "Ljung-Box")
cat("\nLjung-Box test:\n")
cat("  Test statistic:", round(lb_test$statistic, 4), "\n")
cat("  p-value:", round(lb_test$p.value, 4), "\n")
cat("  Interpretation:", 
    ifelse(lb_test$p.value > 0.05, 
           "Residuals are white noise (good)",
           "Residuals show autocorrelation (model may be misspecified)"), "\n")

# 4. Distribution check
cat("\nResidual distribution:\n")
cat("  Skewness:", round(pmm_skewness(residuals_check), 4), "\n")
cat("  Kurtosis:", round(pmm_kurtosis(residuals_check), 4), "\n")

## ----seasonal_adjustment------------------------------------------------------
# Extract seasonal component from fitted model
resid_sar <- as.numeric(fit_sar_pmm2@residuals)
y_fit <- tail(as.numeric(y), length(resid_sar))
fitted_values <- y_fit - resid_sar
seasonal_component <- fitted_values - mean(fitted_values)

# Plot seasonally adjusted series
y_adjusted <- y_fit - seasonal_component

par(mfrow = c(2, 1))
plot(y_fit, type = "l", main = "Original Series (aligned to fitted residuals)",
     col = "steelblue", lwd = 1.5, ylab = "Value")
plot(y_adjusted, type = "l", main = "Seasonally Adjusted Series",
     col = "darkgreen", lwd = 1.5, ylab = "Value")
par(mfrow = c(1, 1))

## ----bootstrap_seasonal, eval=FALSE-------------------------------------------
# # Bootstrap confidence intervals for SAR parameters
# boot_sar <- ts_pmm2_inference(
#   fit_sar_pmm2,
#   B = 500,                    # Number of bootstrap samples
#   method = "block",           # Block bootstrap for time series
#   block_length = 12,          # Block length = seasonal period
#   seed = 2025,
#   parallel = FALSE
# )
# 
# summary(boot_sar)

