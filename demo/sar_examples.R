# SAR-PMM2 Examples and Demonstrations
# ======================================
# Purpose: Demonstrate the usage of SAR (Seasonal AR) models with PMM2 method
# Author: EstemPMM package
# Date: 2025-11-13

library(EstemPMM)

cat("\n")
cat("=", rep("=", 75), "\n", sep = "")
cat("SAR-PMM2: Seasonal Autoregressive Models with Polynomial Maximization\n")
cat("=", rep("=", 75), "\n\n", sep = "")

# ==============================================================================
# EXAMPLE 1: Simple SAR(1)_12 Model
# ==============================================================================

cat("\n")
cat("EXAMPLE 1: Pure Seasonal Model SAR(1)_12\n")
cat("-", rep("-", 75), "\n\n", sep = "")

# Generate synthetic seasonal data with annual pattern
set.seed(123)
n <- 120  # 10 years of monthly data

# Generate SAR(1) with seasonal period 12
# True model: y_t = 0.6 * y_{t-12} + e_t
sar_coef_true <- 0.6
y1 <- numeric(n)

# Initialize first 12 values
y1[1:12] <- rnorm(12, mean = 0, sd = 1)

# Generate rest with seasonal dependence
for (t in 13:n) {
  y1[t] <- sar_coef_true * y1[t - 12] + rnorm(1, sd = 1)
}

# Plot the series
plot(y1, type = "l", main = "Synthetic SAR(1)_12 Series",
     xlab = "Time", ylab = "Value", col = "blue")
abline(h = 0, col = "gray", lty = 2)

# Fit SAR(1)_12 model with PMM2
cat("Fitting SAR(1)_12 model...\n\n")
fit1_pmm2 <- sar_pmm2(y1, order = c(0, 1), season = list(period = 12))

# Display results
summary(fit1_pmm2)

cat("\nTrue SAR coefficient: ", sar_coef_true, "\n")
cat("Estimated SAR coefficient: ", coef(fit1_pmm2)["sar1"], "\n")
cat("Estimation error: ", abs(coef(fit1_pmm2)["sar1"] - sar_coef_true), "\n\n")


# ==============================================================================
# EXAMPLE 2: Combined AR(1) + SAR(1)_12 Model
# ==============================================================================

cat("\n\n")
cat("EXAMPLE 2: Combined AR(1) + SAR(1)_12 Model\n")
cat("-", rep("-", 75), "\n\n", sep = "")

# Generate data with both non-seasonal and seasonal components
set.seed(456)
ar_coef_true <- 0.5
sar_coef_true <- 0.4

y2 <- numeric(n)
y2[1:12] <- rnorm(12, sd = 1)

for (t in 13:n) {
  ar_component <- if (t > 1) ar_coef_true * y2[t - 1] else 0
  sar_component <- sar_coef_true * y2[t - 12]
  y2[t] <- ar_component + sar_component + rnorm(1, sd = 1)
}

# Plot
plot(y2, type = "l", main = "Synthetic AR(1) + SAR(1)_12 Series",
     xlab = "Time", ylab = "Value", col = "darkgreen")
abline(h = 0, col = "gray", lty = 2)

# Fit with PMM2
cat("Fitting AR(1) + SAR(1)_12 model...\n\n")
fit2_pmm2 <- sar_pmm2(y2, order = c(1, 1), season = list(period = 12))

# Display results
summary(fit2_pmm2)

cat("\nTrue coefficients:\n")
cat("  ar1:  ", ar_coef_true, "\n")
cat("  sar1: ", sar_coef_true, "\n")

cat("\nEstimated coefficients:\n")
cat("  ar1:  ", coef(fit2_pmm2)["ar1"], "\n")
cat("  sar1: ", coef(fit2_pmm2)["sar1"], "\n")

cat("\nEstimation errors:\n")
cat("  ar1:  ", abs(coef(fit2_pmm2)["ar1"] - ar_coef_true), "\n")
cat("  sar1: ", abs(coef(fit2_pmm2)["sar1"] - sar_coef_true), "\n\n")


# ==============================================================================
# EXAMPLE 3: SAR with Asymmetric Innovations (Gamma distributed)
# ==============================================================================

cat("\n\n")
cat("EXAMPLE 3: SAR(1)_12 with Asymmetric Innovations\n")
cat("-", rep("-", 75), "\n\n", sep = "")

cat("This example demonstrates PMM2's advantage over OLS when\n")
cat("innovations have asymmetric distribution (right-skewed).\n\n")

# Generate data with gamma-distributed innovations
set.seed(789)
sar_coef_true <- 0.6

y3 <- numeric(n)
y3[1:12] <- rnorm(12, sd = 1)

# Gamma innovations (right-skewed)
gamma_innovations <- rgamma(n, shape = 2, scale = 1) - 2  # Center around 0

for (t in 13:n) {
  y3[t] <- sar_coef_true * y3[t - 12] + gamma_innovations[t]
}

# Plot series and innovation distribution
par(mfrow = c(1, 2))
plot(y3, type = "l", main = "SAR(1)_12 with Gamma Innovations",
     xlab = "Time", ylab = "Value", col = "purple")
hist(gamma_innovations, breaks = 30, main = "Innovation Distribution",
     xlab = "Value", col = "lightblue", border = "white")
par(mfrow = c(1, 1))

# Compare OLS and PMM2
cat("Comparing OLS and PMM2 methods...\n\n")

fit3_ols <- sar_pmm2(y3, order = c(0, 1), season = list(period = 12),
                     method = "ols")
fit3_pmm2 <- sar_pmm2(y3, order = c(0, 1), season = list(period = 12),
                      method = "pmm2")

cat("Results Comparison:\n")
cat("=", rep("=", 75), "\n\n", sep = "")

cat("True SAR coefficient: ", sar_coef_true, "\n\n")

cat("OLS Estimate:\n")
cat("  Coefficient: ", coef(fit3_ols)["sar1"], "\n")
cat("  Error:       ", abs(coef(fit3_ols)["sar1"] - sar_coef_true), "\n")
cat("  Residual SD: ", sqrt(fit3_ols@m2), "\n\n")

cat("PMM2 Estimate:\n")
cat("  Coefficient: ", coef(fit3_pmm2)["sar1"], "\n")
cat("  Error:       ", abs(coef(fit3_pmm2)["sar1"] - sar_coef_true), "\n")
cat("  Residual SD: ", sqrt(fit3_pmm2@m2), "\n\n")

# Calculate improvement
ols_error <- abs(coef(fit3_ols)["sar1"] - sar_coef_true)
pmm2_error <- abs(coef(fit3_pmm2)["sar1"] - sar_coef_true)

if (pmm2_error < ols_error) {
  improvement <- (ols_error - pmm2_error) / ols_error * 100
  cat(sprintf("=> PMM2 reduced estimation error by %.1f%%! [OK]\n\n", improvement))
} else {
  cat("=> OLS performed better in this particular sample\n\n")
}

cat("Distribution characteristics:\n")
cat("  Skewness (c3):       ", fit3_pmm2@m3 / (fit3_pmm2@m2^(3/2)), "\n")
cat("  Excess kurtosis (c4):", fit3_pmm2@m4 / fit3_pmm2@m2^2 - 3, "\n")

c3 <- fit3_pmm2@m3 / (fit3_pmm2@m2^(3/2))
c4 <- fit3_pmm2@m4 / fit3_pmm2@m2^2 - 3
g <- 1 - c3^2 / (2 + c4)

cat("  Variance factor (g): ", g, "\n")
if (g < 1) {
  cat(sprintf("  => Expected %.1f%% variance reduction vs OLS\n\n", (1-g)*100))
}


# ==============================================================================
# EXAMPLE 4: Method Comparison
# ==============================================================================

cat("\n\n")
cat("EXAMPLE 4: Systematic Method Comparison\n")
cat("-", rep("-", 75), "\n\n", sep = "")

cat("Comparing OLS, PMM2, and CSS methods on the same dataset...\n\n")

# Use the asymmetric data from Example 3
compare_sar_methods(y3, order = c(0, 1), period = 12, methods = c("ols", "pmm2", "css"))


# ==============================================================================
# EXAMPLE 5: Quarterly Data (Period = 4)
# ==============================================================================

cat("\n\n")
cat("EXAMPLE 5: Quarterly Data - SAR(1)_4\n")
cat("-", rep("-", 75), "\n\n", sep = "")

# Generate quarterly data (40 quarters = 10 years)
set.seed(321)
n_q <- 40
sar_coef_q <- 0.7

y_quarterly <- numeric(n_q)
y_quarterly[1:4] <- rnorm(4, sd = 1)

for (t in 5:n_q) {
  y_quarterly[t] <- sar_coef_q * y_quarterly[t - 4] + rnorm(1, sd = 0.8)
}

# Plot
plot(y_quarterly, type = "b", pch = 19, col = "darkred",
     main = "Quarterly SAR(1)_4 Series",
     xlab = "Quarter", ylab = "Value")
abline(h = 0, col = "gray", lty = 2)

# Fit SAR(1)_4 model
cat("Fitting SAR(1)_4 model for quarterly data...\n\n")
fit_quarterly <- sar_pmm2(y_quarterly, order = c(0, 1),
                          season = list(period = 4))

summary(fit_quarterly)

cat("\nTrue SAR coefficient: ", sar_coef_q, "\n")
cat("Estimated SAR coefficient: ", coef(fit_quarterly)["sar1"], "\n\n")


# ==============================================================================
# SUMMARY AND CONCLUSIONS
# ==============================================================================

cat("\n\n")
cat("=", rep("=", 75), "\n", sep = "")
cat("SUMMARY AND KEY TAKEAWAYS\n")
cat("=", rep("=", 75), "\n\n", sep = "")

cat("1. SAR-PMM2 successfully models seasonal patterns in time series data\n\n")

cat("2. Key advantages of PMM2 over OLS:\n")
cat("   [OK] More efficient when innovations are asymmetric\n")
cat("   [OK] Lower variance estimates (10-40% reduction possible)\n")
cat("   [OK] Robust to non-Gaussian error distributions\n\n")

cat("3. When to use SAR-PMM2:\n")
cat("   [OK] Monthly data with annual seasonality (period = 12)\n")
cat("   [OK] Quarterly data with annual seasonality (period = 4)\n")
cat("   [OK] Weekly data with annual seasonality (period = 52)\n")
cat("   [OK] Any periodic pattern in the data\n\n")

cat("4. Model selection:\n")
cat("   - Pure seasonal: SAR(0,P)_s\n")
cat("   - With short-term dynamics: SAR(p,P)_s\n")
cat("   - Check ACF/PACF for appropriate orders\n\n")

cat("5. Distribution diagnostics:\n")
cat("   - Check skewness (c3) and kurtosis (c4) of residuals\n")
cat("   - If |c3| > 0.5, expect PMM2 to outperform OLS\n")
cat("   - Variance reduction factor g indicates expected efficiency gain\n\n")

cat("Demo completed successfully!\n\n")

cat("For Monte Carlo simulations demonstrating PMM2 efficiency, run:\n")
cat("  source('demo/sar_monte_carlo.R')\n\n")

cat("For theoretical analysis, see:\n")
cat("  docs/SAR_PMM2_theoretical_analysis.md\n")
cat("  docs/SAR_visual_explanation.md\n\n")
