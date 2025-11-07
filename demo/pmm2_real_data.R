# Demo: PMM2 Applied to Real Data (Auto MPG Dataset)
#
# Purpose: Demonstrate PMM2 on real-world data with actual non-Gaussian residuals
# Duration: 1-2 minutes
#
# This demo shows:
#   - Loading and preprocessing Auto MPG dataset
#   - Comparing OLS vs PMM2 on real data
#   - Residual diagnostics and moment analysis
#   - Visual comparison of fitted models

# Check for required packages
if (!requireNamespace("EstemPMM", quietly = TRUE)) {
  stop("Please install EstemPMM package first", call. = FALSE)
}

library(EstemPMM)

cat("\n")
cat("==============================================================\n")
cat("  PMM2 Applied to Auto MPG Dataset\n")
cat("==============================================================\n\n")

# Load Auto MPG data from UCI Repository
cat("Loading Auto MPG data from UCI Repository...\n")

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/auto-mpg/auto-mpg.data"
auto_data <- try(
  read.table(url, sep = "", strip.white = TRUE, na.strings = "?",
             col.names = c("mpg", "cylinders", "displacement",
                          "horsepower", "weight", "acceleration",
                          "model_year", "origin", "car_name")),
  silent = TRUE
)

if (inherits(auto_data, "try-error")) {
  cat("Could not download data from UCI. Creating synthetic data...\n")
  set.seed(123)
  n <- 200
  x <- runif(n, 8, 25)  # acceleration
  errors <- rgamma(n, shape = 2, scale = 1) - 2  # skewed errors
  y <- 5 + 1.2 * x + errors
  auto_data <- data.frame(mpg = y, acceleration = x)
} else {
  cat("Data loaded successfully!\n")
}

# Preprocessing
auto_data <- auto_data[complete.cases(auto_data[c("mpg", "acceleration")]), ]

cat("Sample size:", nrow(auto_data), "observations\n\n")

# Prepare data
model_data <- data.frame(
  y = auto_data$mpg,
  x = auto_data$acceleration
)

# Fit models
cat("Fitting models...\n")
ols_fit <- lm(y ~ x, data = model_data)
pmm2_fit <- lm_pmm2(y ~ x, data = model_data, verbose = FALSE)

# Extract coefficients
ols_coef <- coef(ols_fit)
pmm2_coef <- coef(pmm2_fit)

# Compute residuals and moments
ols_resid <- residuals(ols_fit)
pmm2_resid <- pmm2_fit@residuals
moments <- compute_moments(ols_resid)

# Bootstrap analysis (OLS vs PMM2 coefficients)
cat("==============================================================\n")
cat("  Bootstrap Analysis (Coefficient Stability)\n")
cat("==============================================================\n\n")

set.seed(2025)
n_boot <- 1000
boot_stat_names <- c("ols_b0", "ols_b1", "pmm2_b0", "pmm2_b1", "diff_b1")
na_boot_vec <- setNames(rep(NA_real_, length(boot_stat_names)), boot_stat_names)
boot_draws <- replicate(
  n_boot,
  {
    idx <- sample.int(nrow(model_data), replace = TRUE)
    dat <- model_data[idx, , drop = FALSE]

    ols_boot <- lm(y ~ x, data = dat)
    pmm2_boot <- try(lm_pmm2(y ~ x, data = dat, verbose = FALSE), silent = TRUE)

    if (inherits(pmm2_boot, "try-error")) {
      return(na_boot_vec)
    }

    ols_coef_b <- coef(ols_boot)
    pmm2_coef_b <- coef(pmm2_boot)

    c(
      ols_b0 = unname(ols_coef_b[1]),
      ols_b1 = unname(ols_coef_b[2]),
      pmm2_b0 = unname(pmm2_coef_b[1]),
      pmm2_b1 = unname(pmm2_coef_b[2]),
      diff_b1 = unname(pmm2_coef_b[2] - ols_coef_b[2])
    )
  },
  simplify = "matrix"
)

rownames(boot_draws) <- boot_stat_names
boot_results <- t(boot_draws)
boot_results <- boot_results[complete.cases(boot_results), , drop = FALSE]

if (nrow(boot_results) == 0) {
  cat("No valid bootstrap samples were obtained. Skipping bootstrap summary.\n\n")
} else {
  cat("Valid resamples:", nrow(boot_results), "of", n_boot, "\n")
  boot_means <- colMeans(boot_results)
  boot_sd <- apply(boot_results, 2, sd)
  boot_ci_b1 <- quantile(boot_results[, "diff_b1"], probs = c(0.025, 0.975))
  boot_ci_b0 <- quantile(
    boot_results[, "pmm2_b0"] - boot_results[, "ols_b0"],
    probs = c(0.025, 0.975)
  )

  cat("Bootstrapped mean slope (beta_1):\n")
  cat("  OLS:  ", sprintf("%.4f", boot_means["ols_b1"]),
      " (SD =", sprintf("%.4f", boot_sd["ols_b1"]), ")\n")
  cat("  PMM2: ", sprintf("%.4f", boot_means["pmm2_b1"]),
      " (SD =", sprintf("%.4f", boot_sd["pmm2_b1"]), ")\n")
  cat("  Difference (PMM2 - OLS):\n")
  cat("    Mean:", sprintf("%+.4f", boot_means["diff_b1"]),
      "| SD =", sprintf("%.4f", boot_sd["diff_b1"]), "\n")
  cat("    95% percentile CI: [",
      sprintf("%+.4f", boot_ci_b1[1]), ", ",
      sprintf("%+.4f", boot_ci_b1[2]), "]\n", sep = "")

  cat("\nBootstrapped intercept difference (PMM2 - OLS):\n")
  cat("  95% percentile CI: [",
      sprintf("%+.4f", boot_ci_b0[1]), ", ",
      sprintf("%+.4f", boot_ci_b0[2]), "]\n\n", sep = "")
}

# Print results
cat("\n==============================================================\n")
cat("  Regression Results\n")
cat("==============================================================\n\n")

cat("OLS Coefficients:\n")
cat("  Intercept:", sprintf("%.4f", ols_coef[1]), "\n")
cat("  Slope:    ", sprintf("%.4f", ols_coef[2]), "\n\n")

cat("PMM2 Coefficients:\n")
cat("  Intercept:", sprintf("%.4f", pmm2_coef[1]), "\n")
cat("  Slope:    ", sprintf("%.4f", pmm2_coef[2]), "\n\n")

cat("Coefficient Differences (PMM2 - OLS):\n")
cat("  Intercept:", sprintf("%+.4f", pmm2_coef[1] - ols_coef[1]), "\n")
cat("  Slope:    ", sprintf("%+.4f", pmm2_coef[2] - ols_coef[2]), "\n\n")

# Model fit statistics
ols_mse <- mean(ols_resid^2)
pmm2_mse <- mean(pmm2_resid^2)
ols_aic <- AIC(ols_fit)
pmm2_aic <- AIC(pmm2_fit)

cat("==============================================================\n")
cat("  Model Fit Statistics\n")
cat("==============================================================\n\n")

cat("Residual MSE:\n")
cat("  OLS:  ", sprintf("%.4f", ols_mse), "\n")
cat("  PMM2: ", sprintf("%.4f", pmm2_mse), "\n")
cat("  Ratio (PMM2/OLS):", sprintf("%.4f", pmm2_mse/ols_mse), "\n\n")

cat("AIC:\n")
cat("  OLS:  ", sprintf("%.2f", ols_aic), "\n")
cat("  PMM2: ", sprintf("%.2f", pmm2_aic), "\n")
cat("  Difference (PMM2 - OLS):", sprintf("%+.2f", pmm2_aic - ols_aic), "\n\n")

# Residual moment analysis
cat("==============================================================\n")
cat("  Residual Distribution Analysis (OLS)\n")
cat("==============================================================\n\n")

cat("Moments:\n")
cat("  Mean (should be ~0):", sprintf("%.4f", mean(ols_resid)), "\n")
cat("  Variance (m2):      ", sprintf("%.4f", moments$m2), "\n")
cat("  Skewness (c3):      ", sprintf("%.4f", moments$c3), "\n")
cat("  Excess Kurtosis (c4):", sprintf("%.4f", moments$c4), "\n\n")

cat("PMM2 Efficiency Indicator:\n")
cat("  g coefficient:", sprintf("%.4f", moments$g), "\n")

if (abs(moments$c3) > 0.5) {
  cat("  -> Residuals show SIGNIFICANT skewness\n")
  cat("  -> PMM2 expected to improve efficiency\n")
} else {
  cat("  -> Residuals close to symmetric\n")
  cat("  -> PMM2 and OLS similar performance expected\n")
}

cat("\n")

# Visualization
cat("Creating diagnostic plots...\n\n")

# Set up 2x2 layout
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. Scatter plot with fitted lines
x_range <- seq(min(model_data$x), max(model_data$x), length.out = 100)
ols_pred <- ols_coef[1] + ols_coef[2] * x_range
pmm2_pred <- pmm2_coef[1] + pmm2_coef[2] * x_range

plot(model_data$x, model_data$y, pch = 16, col = rgb(0, 0, 0, 0.3),
     main = "MPG vs Acceleration",
     xlab = "Acceleration", ylab = "MPG")
lines(x_range, ols_pred, col = "blue", lwd = 2)
lines(x_range, pmm2_pred, col = "red", lwd = 2)
legend("topright", legend = c("OLS", "PMM2", "Data"),
       col = c("blue", "red", rgb(0,0,0,0.3)),
       lty = c(1, 1, NA), pch = c(NA, NA, 16),
       lwd = c(2, 2, NA))

# 2. Histogram of OLS residuals
hist(ols_resid, breaks = 30, probability = TRUE,
     main = "OLS Residuals Distribution",
     xlab = "Residuals", col = "lightblue", border = "white")
lines(density(ols_resid), col = "darkblue", lwd = 2)

# Add normal curve for comparison
x_norm <- seq(min(ols_resid), max(ols_resid), length.out = 100)
y_norm <- dnorm(x_norm, mean = mean(ols_resid), sd = sd(ols_resid))
lines(x_norm, y_norm, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Actual", "Normal"),
       col = c("darkblue", "red"), lty = c(1, 2), lwd = 2, cex = 0.8)

# 3. Q-Q plot
qqnorm(ols_resid, main = "Q-Q Plot (OLS Residuals)", pch = 16, col = rgb(0, 0, 1, 0.5))
qqline(ols_resid, col = "red", lwd = 2)

# 4. Residuals comparison boxplot
boxplot(list(OLS = ols_resid, PMM2 = pmm2_resid),
        main = "Residuals Comparison",
        ylab = "Residual Value",
        col = c("lightblue", "lightgreen"),
        horizontal = FALSE)
abline(h = 0, col = "red", lty = 2)

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Bootstrap visualizations (if available)
if (nrow(boot_results) > 0) {
  cat("Rendering bootstrap coefficient diagnostics...\n\n")
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  hist(
    boot_results[, "diff_b1"],
    breaks = 30,
    probability = TRUE,
    main = expression("Bootstrap Distribution of " * Delta * beta[1]),
    xlab = expression(Delta * beta[1] ~ "(PMM2 - OLS)"),
    col = "lightgreen",
    border = "white"
  )
  lines(density(boot_results[, "diff_b1"]), col = "darkgreen", lwd = 2)
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  abline(v = boot_means["diff_b1"], col = "blue", lwd = 2)
  legend(
    "topright",
    legend = c("Density", "Mean", "Zero difference"),
    col = c("darkgreen", "blue", "red"),
    lty = c(1, 1, 2),
    lwd = c(2, 2, 2),
    cex = 0.8
  )

  plot(ecdf(boot_results[, "diff_b1"]),
       main = expression("ECDF of " * Delta * beta[1]),
       xlab = expression(Delta * beta[1]),
       ylab = "Empirical CDF",
       col = "darkgreen",
       lwd = 2,
       verticals = TRUE,
       do.points = FALSE)
  abline(v = 0, col = "red", lty = 2, lwd = 2)
  abline(v = boot_ci_b1, col = "blue", lty = 3, lwd = 2)
  legend("bottomright",
         legend = c("ECDF", "Zero difference", "95% CI bounds"),
         col = c("darkgreen", "red", "blue"),
         lty = c(1, 2, 3),
         lwd = c(2, 2, 2),
         cex = 0.8)

  boxplot(
    boot_results[, c("ols_b1", "pmm2_b1")],
    names = c("OLS b1", "PMM2 b1"),
    col = c("lightblue", "lightgreen"),
    main = "Bootstrap Slopes",
    ylab = "Slope estimate"
  )
  abline(h = coef(ols_fit)[2], col = "blue", lty = 2)
  abline(h = coef(pmm2_fit)[2], col = "darkgreen", lty = 2)

  boxplot(
    boot_results[, c("ols_b0", "pmm2_b0")],
    names = c("OLS b0", "PMM2 b0"),
    col = c("lightblue", "lightgreen"),
    main = "Bootstrap Intercepts",
    ylab = "Intercept estimate"
  )
  abline(h = coef(ols_fit)[1], col = "blue", lty = 2)
  abline(h = coef(pmm2_fit)[1], col = "darkgreen", lty = 2)

  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

cat("\n==============================================================\n")
cat("  Interpretation\n")
cat("==============================================================\n\n")

cat("1. Coefficient Comparison:\n")
if (abs(pmm2_coef[2] - ols_coef[2])/abs(ols_coef[2]) > 0.05) {
  cat("   - PMM2 and OLS slopes differ by more than 5%\n")
  cat("   - This suggests non-Gaussian errors are present\n")
} else {
  cat("   - PMM2 and OLS slopes are similar\n")
  cat("   - Errors may be approximately Gaussian\n")
}

cat("\n2. Residual Analysis:\n")
cat("   - Q-Q plot shows departures from normality\n")
cat("   - Histogram reveals actual error distribution\n")
if (abs(moments$c3) > 0.5) {
  cat("   - Significant skewness detected (c3 =", sprintf("%.2f", moments$c3), ")\n")
  cat("   - PMM2 should provide more efficient estimates\n")
} else {
  cat("   - Skewness moderate (c3 =", sprintf("%.2f", moments$c3), ")\n")
}

cat("\n3. Model Performance:\n")
if (pmm2_mse < ols_mse) {
  reduction <- (1 - pmm2_mse/ols_mse) * 100
  cat("   - PMM2 achieves", sprintf("%.1f%%", reduction), "lower MSE\n")
  cat("   - Better fit to the data\n")
} else if (pmm2_mse/ols_mse < 1.05) {
  cat("   - PMM2 and OLS have similar MSE\n")
  cat("   - Comparable model fit\n")
} else {
  cat("   - OLS has slightly lower MSE\n")
  cat("   - This can happen with small samples or near-normal errors\n")
}

if (pmm2_aic < ols_aic) {
  cat("   - PMM2 has lower AIC (better model)\n")
} else if (abs(pmm2_aic - ols_aic) < 2) {
  cat("   - AIC values similar (models comparable)\n")
} else {
  cat("   - OLS has lower AIC\n")
}

cat("\n==============================================================\n")
cat("  Conclusion\n")
cat("==============================================================\n\n")

cat("This real-world example demonstrates PMM2 application to\n")
cat("the Auto MPG dataset. The analysis reveals:\n\n")

cat("- PMM2 provides a robust alternative to OLS\n")
cat("- Particularly useful when residual diagnostics show\n")
cat("  departures from normality (skewness, heavy tails)\n")
cat("- Model comparison through MSE and AIC helps select\n")
cat("  the most appropriate method for your data\n\n")

cat("For more detailed statistical inference (confidence intervals,\n")
cat("hypothesis tests), see the package vignettes:\n")
cat("  > browseVignettes(\"EstemPMM\")\n\n")

cat("Demo completed successfully!\n")
cat("==============================================================\n\n")
