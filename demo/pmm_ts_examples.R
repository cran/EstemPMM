# Demo: Time Series Modeling with PMM2
#
# Purpose: Demonstrate PMM2 for AR, MA, ARMA, and ARIMA models
# Duration: 1-2 minutes
#
# This demo shows:
#   - AR(1) and AR(2) models
#   - MA(1) model
#   - ARMA(1,1) model
#   - ARIMA(1,1,1) model
#   - Comparison with classical methods (CSS)
#   - Visual diagnostics

# Check for required package
if (!requireNamespace("EstemPMM", quietly = TRUE)) {
  stop("Please install EstemPMM package first", call. = FALSE)
}

library(EstemPMM)

cat("\n")
cat("==============================================================\n")
cat("  PMM2 for Time Series Models\n")
cat("==============================================================\n\n")

set.seed(2025)

# Helper function to compare methods
compare_ts_fit <- function(series, fit_pmm2, fit_css, model_name, true_coef = NULL) {
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat(" ", model_name, "\n")
  cat("--------------------------------------------------------------\n\n")

  if (!is.null(true_coef)) {
    cat("True coefficients:\n")
    print(true_coef)
    cat("\n")
  }

  cat("CSS estimates:\n")
  print(coef(fit_css))

  cat("\nPMM2 estimates:\n")
  print(coef(fit_pmm2))

  # Residual comparison
  css_resid <- fit_css@residuals
  pmm2_resid <- fit_pmm2@residuals

  cat("\nResidual MSE:\n")
  cat("  CSS:  ", sprintf("%.4f", mean(css_resid^2, na.rm = TRUE)), "\n")
  cat("  PMM2: ", sprintf("%.4f", mean(pmm2_resid^2, na.rm = TRUE)), "\n")

  # Moments
  moments_css <- compute_moments(css_resid[is.finite(css_resid)])
  cat("\nResidual moments (CSS):\n")
  cat("  Skewness (c3):", sprintf("%.4f", moments_css$c3), "\n")
  cat("  Kurtosis (c4):", sprintf("%.4f", moments_css$c4), "\n")

  if (abs(moments_css$c3) > 0.3) {
    cat("  -> Non-Gaussian residuals detected\n")
    cat("  -> PMM2 expected to improve efficiency\n")
  }
}

#===============================================================
# 1. AR(1) MODEL
#===============================================================

cat("\n")
cat("==============================================================\n")
cat("  Example 1: AR(1) Model\n")
cat("==============================================================\n\n")

# Generate AR(1) with skewed innovations
n <- 300
phi1 <- 0.7

# Skewed innovations (gamma distribution)
innovations_ar1 <- rgamma(n, shape = 2, rate = 2) - 1

# Generate AR(1) series
ar1_series <- numeric(n)
ar1_series[1] <- innovations_ar1[1]
for (i in 2:n) {
  ar1_series[i] <- phi1 * ar1_series[i-1] + innovations_ar1[i]
}

cat("Generated AR(1) series with phi =", phi1, "\n")
cat("Sample size:", n, "observations\n")
cat("Innovation distribution: Gamma (right-skewed)\n\n")

# Fit models
cat("Fitting models...\n")
fit_ar1_pmm2 <- ar_pmm2(ar1_series, order = 1, method = "pmm2", include.mean = FALSE)
fit_ar1_css <- ar_pmm2(ar1_series, order = 1, method = "css", include.mean = FALSE)

compare_ts_fit(ar1_series, fit_ar1_pmm2, fit_ar1_css,
               "AR(1) Model", true_coef = phi1)

# Plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Time series plot
plot(ar1_series, type = "l", main = "AR(1) Series",
     xlab = "Time", ylab = "Value", col = "steelblue")

# ACF
acf(ar1_series, main = "ACF of AR(1) Series", lag.max = 20)

# Residuals comparison
boxplot(list(CSS = fit_ar1_css@residuals, PMM2 = fit_ar1_pmm2@residuals),
        main = "Residuals: CSS vs PMM2",
        col = c("lightblue", "lightgreen"))
abline(h = 0, col = "red", lty = 2)

# Q-Q plot of residuals
qqnorm(fit_ar1_pmm2@residuals, main = "Q-Q Plot (PMM2 Residuals)",
       pch = 16, col = rgb(0, 0, 1, 0.5))
qqline(fit_ar1_pmm2@residuals, col = "red", lwd = 2)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

cat("\nPress <Enter> to continue to MA(1) example...")
readline()

#===============================================================
# 2. MA(1) MODEL
#===============================================================

cat("\n")
cat("==============================================================\n")
cat("  Example 2: MA(1) Model\n")
cat("==============================================================\n\n")

# Generate MA(1) with heavy-tailed innovations
theta1 <- 0.6

# Heavy-tailed innovations (t-distribution)
innovations_ma1 <- rt(n, df = 4)

# Generate MA(1) series
ma1_series <- numeric(n)
ma1_series[1] <- innovations_ma1[1]
for (i in 2:n) {
  ma1_series[i] <- innovations_ma1[i] + theta1 * innovations_ma1[i-1]
}

cat("Generated MA(1) series with theta =", theta1, "\n")
cat("Sample size:", n, "observations\n")
cat("Innovation distribution: Student-t (df=4, heavy tails)\n\n")

# Fit models
cat("Fitting models...\n")
fit_ma1_pmm2 <- ma_pmm2(ma1_series, order = 1, method = "pmm2", include.mean = FALSE)
fit_ma1_css <- ma_pmm2(ma1_series, order = 1, method = "css", include.mean = FALSE)

compare_ts_fit(ma1_series, fit_ma1_pmm2, fit_ma1_css,
               "MA(1) Model", true_coef = theta1)

# Plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(ma1_series, type = "l", main = "MA(1) Series",
     xlab = "Time", ylab = "Value", col = "darkgreen")

acf(ma1_series, main = "ACF of MA(1) Series", lag.max = 20)

pacf(ma1_series, main = "PACF of MA(1) Series", lag.max = 20)

# Coefficient comparison
coef_comparison <- rbind(
  True = theta1,
  CSS = coef(fit_ma1_css),
  PMM2 = coef(fit_ma1_pmm2)
)
barplot(coef_comparison, beside = TRUE, main = "MA(1) Coefficient Estimates",
        col = c("gray", "lightblue", "lightgreen"),
        legend.text = TRUE, args.legend = list(x = "topright", cex = 0.8))
abline(h = theta1, col = "red", lty = 2, lwd = 2)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

cat("\nPress <Enter> to continue to ARMA(1,1) example...")
readline()

#===============================================================
# 3. ARMA(1,1) MODEL
#===============================================================

cat("\n")
cat("==============================================================\n")
cat("  Example 3: ARMA(1,1) Model\n")
cat("==============================================================\n\n")

# Generate ARMA(1,1) with mixed innovations
phi_arma <- 0.5
theta_arma <- 0.4

# Mixed distribution: 80% normal + 20% contamination
innovations_arma <- ifelse(runif(n) < 0.8,
                           rnorm(n),
                           rnorm(n, mean = 0, sd = 3))

# Generate ARMA(1,1) series
arma11_series <- arima.sim(n = n,
                           model = list(ar = phi_arma, ma = theta_arma),
                           innov = innovations_arma)

cat("Generated ARMA(1,1) series\n")
cat("  AR coefficient (phi):", phi_arma, "\n")
cat("  MA coefficient (theta):", theta_arma, "\n")
cat("Sample size:", n, "observations\n")
cat("Innovation distribution: Contaminated normal (outliers)\n\n")

# Fit models
cat("Fitting models...\n")
fit_arma_pmm2 <- arma_pmm2(arma11_series, order = c(1, 1),
                           method = "pmm2", include.mean = FALSE)
fit_arma_css <- arma_pmm2(arma11_series, order = c(1, 1),
                          method = "css", include.mean = FALSE)

compare_ts_fit(arma11_series, fit_arma_pmm2, fit_arma_css,
               "ARMA(1,1) Model",
               true_coef = c(ar1 = phi_arma, ma1 = theta_arma))

# Plot
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot(arma11_series, type = "l", main = "ARMA(1,1) Series",
     xlab = "Time", ylab = "Value", col = "purple")

acf(arma11_series, main = "ACF", lag.max = 20)

pacf(arma11_series, main = "PACF", lag.max = 20)

# Residual histogram
hist(fit_arma_pmm2@residuals, breaks = 30, probability = TRUE,
     main = "PMM2 Residuals", xlab = "Residuals",
     col = "lightgreen", border = "white")
lines(density(fit_arma_pmm2@residuals), col = "darkgreen", lwd = 2)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

cat("\nPress <Enter> to continue to ARIMA(1,1,1) example...")
readline()

#===============================================================
# 4. ARIMA(1,1,1) MODEL
#===============================================================

cat("\n")
cat("==============================================================\n")
cat("  Example 4: ARIMA(1,1,1) Model (with differencing)\n")
cat("==============================================================\n\n")

# Generate ARIMA(1,1,1) - integrated series
phi_arima <- 0.4
theta_arima <- 0.3

# Skewed innovations
innovations_arima <- rchisq(n, df = 5) - 5

# Generate stationary ARMA first
arma_stationary <- arima.sim(n = n,
                             model = list(ar = phi_arima, ma = theta_arima),
                             innov = innovations_arima)

# Integrate to create non-stationary series
arima111_series <- cumsum(arma_stationary)

cat("Generated ARIMA(1,1,1) series\n")
cat("  AR coefficient (phi):", phi_arima, "\n")
cat("  MA coefficient (theta):", theta_arima, "\n")
cat("  Differencing order (d): 1\n")
cat("Sample size:", n, "observations\n")
cat("Innovation distribution: Chi-squared (skewed)\n\n")

# Fit models
cat("Fitting models...\n")
fit_arima_pmm2 <- arima_pmm2(arima111_series, order = c(1, 1, 1),
                             method = "pmm2", include.mean = FALSE)
fit_arima_css <- arima_pmm2(arima111_series, order = c(1, 1, 1),
                            method = "css", include.mean = FALSE)

compare_ts_fit(arima111_series, fit_arima_pmm2, fit_arima_css,
               "ARIMA(1,1,1) Model",
               true_coef = c(ar1 = phi_arima, ma1 = theta_arima))

# Plot
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

# Original series (non-stationary)
plot(arima111_series, type = "l", main = "Original Series (Non-stationary)",
     xlab = "Time", ylab = "Value", col = "red")

# Differenced series (stationary)
diff_series <- diff(arima111_series)
plot(diff_series, type = "l", main = "After Differencing (Stationary)",
     xlab = "Time", ylab = "Delta Value", col = "blue")

# ACF of original
acf(arima111_series, main = "ACF (Original)", lag.max = 20)

# ACF of differenced
acf(diff_series, main = "ACF (Differenced)", lag.max = 20)

# PACF of differenced
pacf(diff_series, main = "PACF (Differenced)", lag.max = 20)

# Residuals
boxplot(list(CSS = fit_arima_css@residuals, PMM2 = fit_arima_pmm2@residuals),
        main = "Residuals Comparison",
        col = c("lightblue", "lightgreen"))
abline(h = 0, col = "red", lty = 2)

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

#===============================================================
# SUMMARY
#===============================================================

cat("\n")
cat("==============================================================\n")
cat("  Summary\n")
cat("==============================================================\n\n")

cat("This demo illustrated PMM2 application to various time\n")
cat("series models:\n\n")

cat("1. AR(1) - Autoregressive with skewed innovations\n")
cat("   -> PMM2 provides more efficient estimates\n\n")

cat("2. MA(1) - Moving average with heavy-tailed innovations\n")
cat("   -> PMM2 handles outliers better than CSS\n\n")

cat("3. ARMA(1,1) - Combined model with contaminated innovations\n")
cat("   -> PMM2 robust to mixed distributions\n\n")

cat("4. ARIMA(1,1,1) - Integrated series with skewed innovations\n")
cat("   -> PMM2 works with differenced data\n\n")

cat("Key Observations:\n")
cat("- When innovations are non-Gaussian (skewed, heavy-tailed),\n")
cat("  PMM2 typically provides:\n")
cat("    - Lower residual variance\n")
cat("    - More accurate coefficient estimates\n")
cat("    - Better model fit\n\n")

cat("- For approximately Gaussian innovations, PMM2 and CSS\n")
cat("  perform similarly (as expected)\n\n")

cat("For more advanced analysis:\n")
cat("  - Monte Carlo simulations: demo('pmm2_simMC_ts')\n")
cat("  - Bootstrap inference: vignette('03-bootstrap-inference')\n")
cat("  - Detailed methodology: vignette('02-pmm2-time-series')\n\n")

cat("Demo completed successfully!\n")
cat("==============================================================\n\n")
