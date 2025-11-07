# Quick sanity check for basic EstemPMM functionality
#
# Purpose: Verify that core PMM2 functions work correctly
# Duration: < 30 seconds
#
# This demo shows:
#   - Linear regression: PMM2 vs OLS
#   - AR(1) time series with skewed errors
#   - Basic prediction capabilities

if (!requireNamespace("EstemPMM", quietly = TRUE)) {
  stop("Please install EstemPMM package first", call. = FALSE)
}

library(EstemPMM)

cat("\n")
cat("==============================================================\n")
cat("  Quick PMM2 Sanity Check\n")
cat("==============================================================\n\n")

set.seed(2025)
n <- 150
x1 <- rnorm(n)
x2 <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2
y <- 1.5 + 2 * x1 - 1 * x2 + errors

dat <- data.frame(y = y, x1 = x1, x2 = x2)

cat("--------------------------------------------------------------\n")
cat("  1. Linear Regression: PMM2 vs OLS\n")
cat("--------------------------------------------------------------\n\n")

cat("Fitting models with skewed errors (gamma distribution)...\n\n")

ols_fit <- lm(y ~ x1 + x2, data = dat)
pmm2_fit <- lm_pmm2(y ~ x1 + x2, data = dat)

cat("OLS Summary:\n")
print(summary(ols_fit))

cat("\nPMM2 Summary:\n")
print(summary(pmm2_fit))

cat("\nPMM2 Coefficients:\n")
print(coef(pmm2_fit))

cat("\nTrue coefficients: Intercept = 1.5, x1 = 2.0, x2 = -1.0\n")

# Prediction example
new_data <- data.frame(x1 = c(-1, 0, 1), x2 = c(0.5, -0.3, 1.2))
cat("\nPMM2 Predictions for new data:\n")
print(predict(pmm2_fit, newdata = new_data))

cat("\n")
cat("--------------------------------------------------------------\n")
cat("  2. Time Series: AR(1) with Skewed Errors\n")
cat("--------------------------------------------------------------\n\n")

cat("Generating AR(1) series with gamma-distributed innovations...\n\n")

ts_innov <- rgamma(n, shape = 2, scale = 1) - 2
ts_series <- filter(ts_innov, filter = 0.6, method = "recursive")
ar_fit <- ar_pmm2(ts_series, order = 1)

cat("AR(1) Model Summary:\n")
print(summary(ar_fit))

cat("\nTrue AR coefficient: phi = 0.6\n")

cat("\n")
cat("==============================================================\n")
cat("  Sanity Check Complete\n")
cat("==============================================================\n\n")

cat("ok Linear regression: lm_pmm2() works\n")
cat("ok Time series: ar_pmm2() works\n")
cat("ok Prediction: predict() works\n")
cat("ok Summary: summary() works\n\n")

cat("Next steps:\n")
cat("  - For detailed comparison: demo('pmm2_comparison_boxplots')\n")
cat("  - For real data example: demo('pmm2_real_data')\n")
cat("  - For time series: demo('pmm_ts_examples')\n\n")
