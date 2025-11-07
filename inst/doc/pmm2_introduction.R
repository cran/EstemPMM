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
set.seed(123)

## ----simulate_skewed----------------------------------------------------------
# Sample size
n <- 200

# Generate predictor
x <- rnorm(n, mean = 5, sd = 2)

# True parameters
beta_0 <- 10
beta_1 <- 2.5

# Generate highly skewed errors using chi-squared distribution
# (chi-square with df=3 has skewness ~ 1.63)
errors <- rchisq(n, df = 3) - 3  # Center to have mean ~ 0

# Response variable
y <- beta_0 + beta_1 * x + errors

# Create data frame
dat <- data.frame(y = y, x = x)

## ----plot_errors, fig.width=7, fig.height=4-----------------------------------
# Fit OLS to get residuals for visualization
fit_ols_temp <- lm(y ~ x, data = dat)
residuals_ols <- residuals(fit_ols_temp)

# Plot histogram of residuals
par(mfrow = c(1, 2))
hist(residuals_ols, breaks = 30, main = "Residual Distribution",
     xlab = "Residuals", col = "lightblue", border = "white")
qqnorm(residuals_ols, main = "Q-Q Plot")
qqline(residuals_ols, col = "red", lwd = 2)
par(mfrow = c(1, 1))

# Calculate skewness and kurtosis
cat(sprintf("Skewness: %.3f\n", pmm_skewness(residuals_ols)))
cat(sprintf("Excess Kurtosis: %.3f\n", pmm_kurtosis(residuals_ols) - 3))

## ----fit_models---------------------------------------------------------------
# Fit PMM2 model
fit_pmm2 <- lm_pmm2(y ~ x, data = dat, verbose = FALSE)

# Fit OLS model for comparison
fit_ols <- lm(y ~ x, data = dat)

# Display coefficients
cat("=== OLS Estimates ===\n")
print(coef(fit_ols))

cat("\n=== PMM2 Estimates ===\n")
print(coef(fit_pmm2))

cat("\n=== True Values ===\n")
cat("Intercept:", beta_0, "\n")
cat("Slope:    ", beta_1, "\n")

## ----bootstrap_inference------------------------------------------------------
# Perform bootstrap inference
boot_results <- pmm2_inference(fit_pmm2, formula = y ~ x, data = dat,
                               B = 500)

# Display summary with confidence intervals
summary(boot_results)

## ----multiple_regression------------------------------------------------------
# Generate multiple predictors
n <- 250
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- runif(n, min = -2, max = 2)
x3 <- rexp(n, rate = 0.5)

# True model
beta <- c(5, 1.2, -0.8, 2.0)  # Intercept, x1, x2, x3
errors <- rt(n, df = 4)  # t-distribution with heavy tails

y <- beta[1] + beta[2]*x1 + beta[3]*x2 + beta[4]*x3 + errors

# Create data frame
dat_multi <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

# Fit models
fit_pmm2_multi <- lm_pmm2(y ~ x1 + x2 + x3, data = dat_multi, verbose = FALSE)
fit_ols_multi <- lm(y ~ x1 + x2 + x3, data = dat_multi)

# Compare coefficients
comparison <- data.frame(
  Parameter = c("Intercept", "x1", "x2", "x3"),
  True = beta,
  OLS = coef(fit_ols_multi),
  PMM2 = coef(fit_pmm2_multi)
)

print(comparison, row.names = FALSE)

## ----polynomial_regression----------------------------------------------------
# Generate data with quadratic relationship
n <- 200
x <- seq(-3, 3, length.out = n)
errors <- rchisq(n, df = 5) - 5

# True quadratic model: y = 2 + 3x - 0.5x^2
y <- 2 + 3*x - 0.5*x^2 + errors

dat_poly <- data.frame(y = y, x = x)

# Fit polynomial model
fit_pmm2_poly <- lm_pmm2(y ~ x + I(x^2), data = dat_poly, verbose = FALSE)

# Summary
summary(fit_pmm2_poly, formula = y ~ x + I(x^2), data = dat_poly)

## ----diagnostics, fig.width=7, fig.height=6-----------------------------------
# Generate diagnostic plots
plot(fit_pmm2_multi)
title("PMM2 Diagnostics")

