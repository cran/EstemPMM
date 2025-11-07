# Demo: Prediction Accuracy Comparison (PMM2 vs OLS)
#
# Purpose: Compare out-of-sample prediction performance using train/test split
# Duration: 1-2 minutes
#
# This demo shows:
#   - Train/test data splitting
#   - Model training on training set
#   - Prediction on unseen test data
#   - Comparison of prediction errors (MSE, MAE, R^2)

# Check for required package
if (!requireNamespace("EstemPMM", quietly = TRUE)) {
  stop("Please install EstemPMM package first", call. = FALSE)
}

library(EstemPMM)

cat("\n")
cat("==============================================================\n")
cat("  PMM2 vs OLS: Prediction Accuracy Comparison\n")
cat("==============================================================\n\n")

# Set seed for reproducibility
set.seed(42)

# Generate synthetic data with skewed errors
cat("Generating synthetic data with skewed errors...\n")
n <- 300
x <- rnorm(n, mean = 5, sd = 2)

# True parameters
true_beta0 <- 10
true_beta1 <- 2.5

# Generate skewed errors (chi-squared distribution)
errors <- rchisq(n, df = 4) - 4

# Response variable
y <- true_beta0 + true_beta1 * x + errors

# Create data frame
full_data <- data.frame(x = x, y = y)

cat("Total sample size:", n, "observations\n")
cat("True parameters: intercept =", true_beta0, ", slope =", true_beta1, "\n")
cat("Error distribution: chi^2(df=4) - 4 (right-skewed)\n\n")

# Split data into training (80%) and test (20%) sets
cat("Splitting data: 80% training, 20% test...\n")
train_size <- floor(0.8 * n)
train_indices <- sample(1:n, train_size)

train_data <- full_data[train_indices, ]
test_data <- full_data[-train_indices, ]

cat("Training set:", nrow(train_data), "observations\n")
cat("Test set:    ", nrow(test_data), "observations\n\n")

# Train models on training data
cat("Training models on training data...\n")
ols_fit <- lm(y ~ x, data = train_data)
pmm2_fit <- lm_pmm2(y ~ x, data = train_data, verbose = FALSE)

cat("Training completed!\n\n")

# Print training set coefficients
cat("==============================================================\n")
cat("  Coefficients (Estimated on Training Data)\n")
cat("==============================================================\n\n")

ols_coef <- coef(ols_fit)
pmm2_coef <- coef(pmm2_fit)

cat("True values:\n")
cat("  Intercept:", sprintf("%.4f", true_beta0), "\n")
cat("  Slope:    ", sprintf("%.4f", true_beta1), "\n\n")

cat("OLS estimates:\n")
cat("  Intercept:", sprintf("%.4f", ols_coef[1]),
    sprintf("(error: %+.4f)", ols_coef[1] - true_beta0), "\n")
cat("  Slope:    ", sprintf("%.4f", ols_coef[2]),
    sprintf("(error: %+.4f)", ols_coef[2] - true_beta1), "\n\n")

cat("PMM2 estimates:\n")
cat("  Intercept:", sprintf("%.4f", pmm2_coef[1]),
    sprintf("(error: %+.4f)", pmm2_coef[1] - true_beta0), "\n")
cat("  Slope:    ", sprintf("%.4f", pmm2_coef[2]),
    sprintf("(error: %+.4f)", pmm2_coef[2] - true_beta1), "\n\n")

# Make predictions on test data
cat("==============================================================\n")
cat("  Prediction Performance (Test Data)\n")
cat("==============================================================\n\n")

cat("Making predictions on test set...\n")
ols_pred <- predict(ols_fit, newdata = test_data)
pmm2_pred <- predict(pmm2_fit, newdata = test_data)

# Calculate prediction errors
ols_errors <- test_data$y - ols_pred
pmm2_errors <- test_data$y - pmm2_pred

# Calculate performance metrics
ols_mse <- mean(ols_errors^2)
pmm2_mse <- mean(pmm2_errors^2)

ols_mae <- mean(abs(ols_errors))
pmm2_mae <- mean(abs(pmm2_errors))

# R-squared on test data
ss_total <- sum((test_data$y - mean(test_data$y))^2)
ols_r2 <- 1 - sum(ols_errors^2) / ss_total
pmm2_r2 <- 1 - sum(pmm2_errors^2) / ss_total

# Print metrics
cat("\nMean Squared Error (MSE):\n")
cat("  OLS:  ", sprintf("%.4f", ols_mse), "\n")
cat("  PMM2: ", sprintf("%.4f", pmm2_mse), "\n")
cat("  Improvement:", sprintf("%.2f%%", (1 - pmm2_mse/ols_mse) * 100), "\n\n")

cat("Mean Absolute Error (MAE):\n")
cat("  OLS:  ", sprintf("%.4f", ols_mae), "\n")
cat("  PMM2: ", sprintf("%.4f", pmm2_mae), "\n")
cat("  Improvement:", sprintf("%.2f%%", (1 - pmm2_mae/ols_mae) * 100), "\n\n")

cat("R-squared (on test data):\n")
cat("  OLS:  ", sprintf("%.4f", ols_r2), "\n")
cat("  PMM2: ", sprintf("%.4f", pmm2_r2), "\n\n")

# Visualization
cat("Creating visualization...\n\n")

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. Training data with fitted lines
plot(train_data$x, train_data$y, pch = 16, col = rgb(0, 0.5, 0, 0.3),
     main = "Training Data with Fitted Models",
     xlab = "x", ylab = "y")
x_range <- seq(min(full_data$x), max(full_data$x), length.out = 100)
lines(x_range, ols_coef[1] + ols_coef[2] * x_range, col = "blue", lwd = 2)
lines(x_range, pmm2_coef[1] + pmm2_coef[2] * x_range, col = "red", lwd = 2)
abline(a = true_beta0, b = true_beta1, col = "black", lty = 2, lwd = 2)
legend("topleft", legend = c("OLS", "PMM2", "True", "Train data"),
       col = c("blue", "red", "black", rgb(0,0.5,0,0.3)),
       lty = c(1, 1, 2, NA), pch = c(NA, NA, NA, 16),
       lwd = c(2, 2, 2, NA), cex = 0.8)

# 2. Test data with predictions
plot(test_data$x, test_data$y, pch = 16, col = rgb(1, 0, 0, 0.5),
     main = "Test Data with Predictions",
     xlab = "x", ylab = "y")
lines(x_range, ols_coef[1] + ols_coef[2] * x_range, col = "blue", lwd = 2)
lines(x_range, pmm2_coef[1] + pmm2_coef[2] * x_range, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("OLS", "PMM2", "Test data"),
       col = c("blue", "red", rgb(1,0,0,0.5)),
       lty = c(1, 2, NA), pch = c(NA, NA, 16),
       lwd = c(2, 2, NA), cex = 0.8)

# 3. Prediction errors comparison (boxplot)
boxplot(list(OLS = ols_errors, PMM2 = pmm2_errors),
        main = "Prediction Errors on Test Data",
        ylab = "Prediction Error",
        col = c("lightblue", "lightgreen"))
abline(h = 0, col = "red", lty = 2)

# 4. Predicted vs Actual
plot_range <- range(c(test_data$y, ols_pred, pmm2_pred))
plot(test_data$y, ols_pred, pch = 16, col = rgb(0, 0, 1, 0.5),
     main = "Predicted vs Actual (Test Data)",
     xlab = "Actual", ylab = "Predicted",
     xlim = plot_range, ylim = plot_range)
points(test_data$y, pmm2_pred, pch = 17, col = rgb(1, 0, 0, 0.5))
abline(a = 0, b = 1, col = "black", lty = 2)
legend("topleft", legend = c("OLS", "PMM2", "Perfect"),
       col = c(rgb(0,0,1,0.5), rgb(1,0,0,0.5), "black"),
       pch = c(16, 17, NA), lty = c(NA, NA, 2), cex = 0.8)

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

# Summary table
cat("==============================================================\n")
cat("  Summary Table\n")
cat("==============================================================\n\n")

summary_df <- data.frame(
  Metric = c("MSE", "MAE", "R^2", "MSE Ratio"),
  OLS = c(ols_mse, ols_mae, ols_r2, 1.0),
  PMM2 = c(pmm2_mse, pmm2_mae, pmm2_r2, pmm2_mse/ols_mse),
  Improvement = c(
    sprintf("%.2f%%", (1 - pmm2_mse/ols_mse) * 100),
    sprintf("%.2f%%", (1 - pmm2_mae/ols_mae) * 100),
    sprintf("%+.4f", pmm2_r2 - ols_r2),
    sprintf("%.4f", pmm2_mse/ols_mse)
  )
)

print(summary_df, row.names = FALSE)

cat("\n==============================================================\n")
cat("  Interpretation\n")
cat("==============================================================\n\n")

if (pmm2_mse < ols_mse) {
  improvement <- (1 - pmm2_mse/ols_mse) * 100
  cat("PMM2 outperforms OLS on out-of-sample prediction:\n")
  cat("  - MSE reduced by", sprintf("%.2f%%", improvement), "\n")
  cat("  - This demonstrates PMM2's robustness to non-Gaussian errors\n")
  cat("  - Better generalization to unseen data\n")
} else {
  cat("OLS and PMM2 show comparable prediction performance:\n")
  cat("  - This can occur when:\n")
  cat("    - Test set happens to have different error characteristics\n")
  cat("    - Sample size is small\n")
  cat("    - Errors are closer to Gaussian in test set\n")
}

cat("\nKey Observations:\n")
cat("  - Both methods trained on the same training data\n")
cat("  - Predictions evaluated on completely unseen test data\n")
cat("  - PMM2 leverages higher-order moments for better estimation\n")
cat("  - With skewed errors, PMM2 typically provides better\n")
cat("    out-of-sample predictions\n")

cat("\n==============================================================\n")
cat("  Conclusion\n")
cat("==============================================================\n\n")

cat("This demo illustrates the importance of testing model\n")
cat("performance on held-out data. PMM2's efficiency gains\n")
cat("translate to:\n\n")

cat("  - Lower prediction errors on new data\n")
cat("  - Better model generalization\n")
cat("  - Improved reliability in production settings\n\n")

cat("For production applications with non-Gaussian errors,\n")
cat("PMM2 offers a principled approach to achieving better\n")
cat("predictive performance.\n\n")

cat("Demo completed successfully!\n")
cat("==============================================================\n\n")
