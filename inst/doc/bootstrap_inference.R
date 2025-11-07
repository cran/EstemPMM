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
set.seed(2025)

## ----basic_linear_bootstrap---------------------------------------------------
# Generate data with skewed errors
n <- 200
x <- rnorm(n, mean = 5, sd = 2)

# True parameters
beta_0 <- 10
beta_1 <- 2.5

# Skewed errors: chi-squared distribution
errors <- rchisq(n, df = 4) - 4

# Response
y <- beta_0 + beta_1 * x + errors
dat <- data.frame(y = y, x = x)

# Fit PMM2 model
fit_pmm2 <- lm_pmm2(y ~ x, data = dat, verbose = FALSE)

# Perform bootstrap inference
boot_results <- pmm2_inference(fit_pmm2,
                               formula = y ~ x,
                               data = dat,
                               B = 1000,
                               seed = 123)

# Display summary
summary(boot_results)

## ----bootstrap_distribution, fig.width=7, fig.height=6------------------------
# Extract bootstrap samples for visualization
# (Note: This requires accessing internals - in practice, use summary output)

# Refit to get bootstrap samples
set.seed(123)
B <- 1000
boot_samples <- matrix(0, nrow = B, ncol = 2)

for(b in 1:B) {
  # Resample residuals
  res_b <- sample(residuals(fit_pmm2), size = n, replace = TRUE)

  # Generate bootstrap data
  y_b <- fitted(fit_pmm2) + res_b
  dat_b <- dat
  dat_b$y <- y_b

  # Refit
  fit_b <- lm_pmm2(y ~ x, data = dat_b, verbose = FALSE)
  boot_samples[b, ] <- coef(fit_b)
}

# Plot bootstrap distributions
par(mfrow = c(2, 2))

# Intercept
hist(boot_samples[, 1], breaks = 30, main = "Bootstrap: Intercept",
     xlab = "Intercept", col = "lightblue", border = "white")
abline(v = beta_0, col = "red", lwd = 2, lty = 2)
abline(v = coef(fit_pmm2)[1], col = "darkblue", lwd = 2)
legend("topright", legend = c("True", "Estimate"),
       col = c("red", "darkblue"), lty = c(2, 1), lwd = 2, cex = 0.8)

# Q-Q plot for intercept
qqnorm(boot_samples[, 1], main = "Q-Q Plot: Intercept")
qqline(boot_samples[, 1], col = "red", lwd = 2)

# Slope
hist(boot_samples[, 2], breaks = 30, main = "Bootstrap: Slope",
     xlab = "Slope", col = "lightgreen", border = "white")
abline(v = beta_1, col = "red", lwd = 2, lty = 2)
abline(v = coef(fit_pmm2)[2], col = "darkblue", lwd = 2)
legend("topright", legend = c("True", "Estimate"),
       col = c("red", "darkblue"), lty = c(2, 1), lwd = 2, cex = 0.8)

# Q-Q plot for slope
qqnorm(boot_samples[, 2], main = "Q-Q Plot: Slope")
qqline(boot_samples[, 2], col = "red", lwd = 2)

par(mfrow = c(1, 1))

## ----ci_comparison------------------------------------------------------------
# Fit OLS for comparison
fit_ols <- lm(y ~ x, data = dat)

# Extract standard errors and CIs
se_ols <- summary(fit_ols)$coefficients[, "Std. Error"]
ci_ols <- confint(fit_ols)

# Bootstrap CIs (from summary)
boot_summary <- boot_results

# Create comparison table
comparison <- data.frame(
  Parameter = c("Intercept", "Slope"),
  True = c(beta_0, beta_1),
  PMM2_Estimate = coef(fit_pmm2),
  OLS_SE = se_ols,
  PMM2_Boot_SE = boot_summary$Std.Error,
  SE_Ratio = boot_summary$Std.Error / se_ols
)

print(comparison, row.names = FALSE, digits = 3)

cat("\n=== Confidence Interval Comparison ===\n\n")
cat("OLS 95% CI (Normal-based):\n")
print(ci_ols)

cat("\nPMM2 95% CI (Bootstrap Percentile):\n")
print(data.frame(
  Parameter = rownames(ci_ols),
  Lower = boot_summary$conf.low,
  Upper = boot_summary$conf.high
))

## ----hypothesis_testing-------------------------------------------------------
# H0: Slope = 0 (no relationship)
# From bootstrap summary
boot_summary <- boot_results

cat("=== Hypothesis Test: Slope = 0 ===\n\n")
cat("t-statistic:", boot_summary$t.value[2], "\n")
cat("p-value:    ", boot_summary$p.value[2], "\n\n")

if(boot_summary$p.value[2] < 0.05) {
  cat("Decision: Reject H0 at 5% level\n")
  cat("Conclusion: Significant relationship between x and y\n")
} else {
  cat("Decision: Fail to reject H0\n")
}

## ----multiple_regression_test-------------------------------------------------
# Generate multiple regression data
n <- 250
x1 <- rnorm(n)
x2 <- rnorm(n)
errors <- rexp(n, rate = 1) - 1

# True model: similar coefficients
beta <- c(5, 1.5, 1.8, -0.5)
y <- beta[1] + beta[2]*x1 + beta[3]*x2 + beta[4]*x1*x2 + errors

dat_multi <- data.frame(y = y, x1 = x1, x2 = x2)

# Fit model
fit_multi <- lm_pmm2(y ~ x1 + x2 + x1:x2, data = dat_multi, verbose = FALSE)

# Bootstrap
boot_multi <- pmm2_inference(fit_multi,
                             formula = y ~ x1 + x2 + x1:x2,
                             data = dat_multi,
                             B = 1000,
                             seed = 456)

# Test H0: beta_x1 = beta_x2
summary(boot_multi)

## ----ar_bootstrap-------------------------------------------------------------
# Generate AR(1) data with skewed innovations
n <- 300
phi <- 0.7
innovations <- rgamma(n, shape = 2, rate = 2) - 1

x <- numeric(n)
x[1] <- innovations[1]
for(i in 2:n) {
  x[i] <- phi * x[i-1] + innovations[i]
}

# Fit AR(1) model
fit_ar <- ar_pmm2(x, order = 1, method = "pmm2", include.mean = FALSE)

# Bootstrap inference for time series
boot_ar <- ts_pmm2_inference(fit_ar,
                             B = 1000,
                             seed = 789)

# Summary
summary(boot_ar)

# Check if AR coefficient is significantly different from 0.5
cat("\n=== Testing H0: phi = 0.5 ===\n")
boot_summary_ar <- boot_ar
ci_lower <- boot_summary_ar$conf.low[1]
ci_upper <- boot_summary_ar$conf.high[1]

cat("95% CI: [", round(ci_lower, 3), ",", round(ci_upper, 3), "]\n")

if(0.5 >= ci_lower && 0.5 <= ci_upper) {
  cat("Decision: Cannot reject H0 (0.5 is within CI)\n")
} else {
  cat("Decision: Reject H0 (0.5 is outside CI)\n")
}

## ----arma_bootstrap-----------------------------------------------------------
# Generate ARMA(1,1) with heavy-tailed innovations
n <- 400
phi <- 0.5
theta <- 0.4
innovations <- rt(n, df = 4)

x <- arima.sim(n = n,
               model = list(ar = phi, ma = theta),
               innov = innovations)

# Fit ARMA model
fit_arma <- arma_pmm2(x, order = c(1, 1), method = "pmm2", include.mean = FALSE)

# Bootstrap
boot_arma <- ts_pmm2_inference(fit_arma,
                               B = 1000,
                               seed = 999,
                               method = "block",
                               block_length = 20)

# Display results
summary(boot_arma)

## ----bootstrap_stability, eval=FALSE------------------------------------------
# # Assess stability of bootstrap SEs with different B
# B_values <- c(100, 500, 1000, 2000)
# se_results <- data.frame(B = B_values, SE_Intercept = NA, SE_Slope = NA)
# 
# for(i in seq_along(B_values)) {
#   boot_temp <- pmm2_inference(fit_pmm2, formula = y ~ x, data = dat,
#                               B = B_values[i], seed = 123)
#   summary_temp <- boot_temp
#   se_results$SE_Intercept[i] <- summary_temp$Std.Error[1]
#   se_results$SE_Slope[i] <- summary_temp$Std.Error[2]
# }
# 
# print(se_results)

## ----parallel_bootstrap, eval=FALSE-------------------------------------------
# # Enable parallel bootstrap (requires 'parallel' package)
# boot_parallel <- pmm2_inference(fit_pmm2,
#                                 formula = y ~ x,
#                                 data = dat,
#                                 B = 2000,
#                                 seed = 123,
#                                 parallel = TRUE,
#                                 cores = 4)  # Use 4 CPU cores
# 
# # Check speedup
# # Typical speedup: 3-4x on 4 cores

## ----bootstrap_diagnostics----------------------------------------------------
# Rerun bootstrap with verbose output for diagnostics
boot_verbose <- pmm2_inference(fit_pmm2,
                               formula = y ~ x,
                               data = dat,
                               B = 500,
                               seed = 321)

# Check for failed bootstrap samples
# (In production code, inspect convergence attribute)
summary(boot_verbose)

cat("\nAll bootstrap samples should converge successfully.\n")
cat("If many fail, consider:\n")
cat("- Increasing max_iter in lm_pmm2()\n")
cat("- Checking for outliers or influential points\n")
cat("- Simplifying the model\n")

## ----boot_viz, fig.width=7, fig.height=4--------------------------------------
# For demonstration, regenerate bootstrap samples
set.seed(123)
B <- 1000
sample_size <- nrow(dat)
boot_coefs <- matrix(0, nrow = B, ncol = 2)

for(b in 1:B) {
  res_b <- sample(residuals(fit_pmm2), sample_size, replace = TRUE)
  y_b <- fitted(fit_pmm2) + res_b
  dat_b <- dat
  dat_b$y <- y_b

  fit_b <- lm_pmm2(y ~ x, data = dat_b, verbose = FALSE)
  boot_coefs[b, ] <- coef(fit_b)
}

# Scatter plot of bootstrap estimates
par(mfrow = c(1, 2))

# Bivariate distribution
plot(boot_coefs[, 1], boot_coefs[, 2],
     pch = 16, col = rgb(0, 0, 1, 0.1),
     xlab = "Intercept", ylab = "Slope",
     main = "Joint Bootstrap Distribution")
points(coef(fit_pmm2)[1], coef(fit_pmm2)[2],
       pch = 19, col = "red", cex = 2)

# Correlation
cor_boot <- cor(boot_coefs[, 1], boot_coefs[, 2])
cat("Bootstrap correlation between intercept and slope:", round(cor_boot, 3), "\n")

# Density plot
plot(density(boot_coefs[, 2]), main = "Slope Bootstrap Density",
     xlab = "Slope", lwd = 2, col = "blue")
abline(v = beta_1, col = "red", lwd = 2, lty = 2)
abline(v = coef(fit_pmm2)[2], col = "darkblue", lwd = 2)
legend("topright", legend = c("True", "Estimate", "Bootstrap"),
       col = c("red", "darkblue", "blue"),
       lty = c(2, 1, 1), lwd = 2)

par(mfrow = c(1, 1))

## ----studentized_bootstrap, eval=FALSE----------------------------------------
# # Not yet implemented in EstemPMM
# # Future feature: studentized CIs for improved small-sample performance

## ----block_bootstrap, eval=FALSE----------------------------------------------
# # Not yet implemented in EstemPMM
# # Current approach: residual bootstrap (adequate for most cases)
# # Future: moving block bootstrap for highly dependent series

