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

## ----simulate_uniform---------------------------------------------------------
n <- 300

# Generate predictor
x <- rnorm(n, mean = 5, sd = 2)

# True parameters
beta_0 <- 10
beta_1 <- 2.5

# Generate uniform errors (platykurtic, gamma4 ≈ -1.2)
errors <- runif(n, -sqrt(3), sqrt(3))  # variance = 1

# Response variable
y <- beta_0 + beta_1 * x + errors
dat <- data.frame(y = y, x = x)

## ----plot_errors_uniform, fig.width=7, fig.height=4---------------------------
fit_ols_temp <- lm(y ~ x, data = dat)
res_ols <- residuals(fit_ols_temp)

par(mfrow = c(1, 2))
hist(res_ols, breaks = 30, main = "Residual Distribution (Uniform Errors)",
     xlab = "Residuals", col = "lightblue", border = "white")
qqnorm(res_ols, main = "Q-Q Plot")
qqline(res_ols, col = "red", lwd = 2)
par(mfrow = c(1, 1))

cat("Skewness (gamma3):", round(pmm_skewness(res_ols), 3), "\n")
cat("Excess kurtosis (gamma4):", round(pmm_kurtosis(res_ols) - 3, 3), "\n")

## ----dispatch_uniform---------------------------------------------------------
recommendation <- pmm_dispatch(res_ols)

## ----fit_pmm3_simple----------------------------------------------------------
fit_pmm3 <- lm_pmm3(y ~ x, data = dat)
fit_ols  <- lm(y ~ x, data = dat)

# Display PMM3 summary
summary(fit_pmm3)

## ----compare_simple-----------------------------------------------------------
comparison <- data.frame(
  Parameter = c("Intercept", "Slope"),
  True  = c(beta_0, beta_1),
  OLS   = coef(fit_ols),
  PMM3  = coef(fit_pmm3),
  Diff  = coef(fit_pmm3) - coef(fit_ols)
)
print(comparison, row.names = FALSE)

cat("\nResidual sum of squares:\n")
cat("  OLS:  ", sum(residuals(fit_ols)^2), "\n")
cat("  PMM3: ", sum(residuals(fit_pmm3)^2), "\n")

cat("\nTheoretical variance ratio g3 =",
    fit_pmm3@g_coefficient, "\n")
cat("Expected variance reduction:",
    round((1 - fit_pmm3@g_coefficient) * 100, 1), "%\n")

## ----three_way----------------------------------------------------------------
# Fit all three methods
fit_ols  <- lm(y ~ x, data = dat)
fit_pmm2 <- lm_pmm2(y ~ x, data = dat)
fit_pmm3 <- lm_pmm3(y ~ x, data = dat)

comparison3 <- data.frame(
  Method = c("OLS", "PMM2", "PMM3"),
  Intercept = c(coef(fit_ols)[1], coef(fit_pmm2)[1], coef(fit_pmm3)[1]),
  Slope     = c(coef(fit_ols)[2], coef(fit_pmm2)[2], coef(fit_pmm3)[2])
)
print(comparison3, row.names = FALSE, digits = 5)

cat("\nTrue values: Intercept =", beta_0, ", Slope =", beta_1, "\n")
cat("\nEfficiency factors:\n")
vf2 <- pmm2_variance_factor(fit_pmm2@m2, fit_pmm2@m3, fit_pmm2@m4)
cat("  PMM2 g2 =", vf2$g2,
    " (improvement:", round((1 - vf2$g2) * 100, 1), "%)\n")
cat("  PMM3 g3 =", fit_pmm3@g_coefficient,
    " (improvement:", round((1 - fit_pmm3@g_coefficient) * 100, 1), "%)\n")

## ----auto_mpg_load------------------------------------------------------------
data(auto_mpg)

# Remove missing horsepower values
auto_complete <- na.omit(auto_mpg[, c("mpg", "horsepower")])

# Visualize the relationship (clearly nonlinear)
plot(auto_complete$horsepower, auto_complete$mpg,
     main = "MPG vs Horsepower",
     xlab = "Horsepower", ylab = "Miles per Gallon",
     col = "steelblue", pch = 16, cex = 0.8)

## ----auto_mpg_fit-------------------------------------------------------------
# Quadratic OLS
fit_auto_ols <- lm(mpg ~ horsepower + I(horsepower^2), data = auto_complete)

# Check residual properties
res_auto <- residuals(fit_auto_ols)
cat("Residual diagnostics:\n")
cat("  Skewness (gamma3):", round(pmm_skewness(res_auto), 3), "\n")
cat("  Kurtosis (gamma4):", round(pmm_kurtosis(res_auto) - 3, 3), "\n")
sym_test <- test_symmetry(res_auto)
cat("  Symmetric:", sym_test$is_symmetric, "\n")

# Dispatch recommendation
dispatch_auto <- pmm_dispatch(res_auto)

## ----auto_mpg_pmm3------------------------------------------------------------
# Fit PMM3 (quadratic model)
fit_auto_pmm3 <- lm_pmm3(mpg ~ horsepower + I(horsepower^2),
                          data = auto_complete)

# Fit PMM2 for comparison
fit_auto_pmm2 <- lm_pmm2(mpg ~ horsepower + I(horsepower^2),
                          data = auto_complete, na.action = na.omit)

# Compare all three methods
comparison_auto <- data.frame(
  Method    = c("OLS", "PMM2", "PMM3"),
  Intercept = c(coef(fit_auto_ols)[1], coef(fit_auto_pmm2)[1],
                coef(fit_auto_pmm3)[1]),
  HP        = c(coef(fit_auto_ols)[2], coef(fit_auto_pmm2)[2],
                coef(fit_auto_pmm3)[2]),
  HP_sq     = c(coef(fit_auto_ols)[3], coef(fit_auto_pmm2)[3],
                coef(fit_auto_pmm3)[3])
)
print(comparison_auto, row.names = FALSE, digits = 6)

cat("\nPMM3 g3 =", fit_auto_pmm3@g_coefficient,
    " (improvement:", round((1 - fit_auto_pmm3@g_coefficient) * 100, 1), "%)\n")

## ----auto_mpg_plot, fig.width=7, fig.height=5---------------------------------
hp_seq <- seq(min(auto_complete$horsepower),
              max(auto_complete$horsepower), length.out = 200)
newdata <- data.frame(horsepower = hp_seq)

pred_ols  <- predict(fit_auto_ols, newdata = newdata)
pred_pmm3 <- predict(fit_auto_pmm3,
                     newdata = data.frame(horsepower = hp_seq,
                                          `I(horsepower^2)` = hp_seq^2))

plot(auto_complete$horsepower, auto_complete$mpg,
     main = "MPG vs Horsepower: OLS and PMM3 Quadratic Fits",
     xlab = "Horsepower", ylab = "Miles per Gallon",
     col = "gray70", pch = 16, cex = 0.7)
lines(hp_seq, pred_ols, col = "blue", lwd = 2)
lines(hp_seq, as.numeric(pred_pmm3), col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("OLS", "PMM3"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

## ----auto_mpg_acceleration----------------------------------------------------
# Linear model: MPG vs Acceleration
fit_accel_ols <- lm(mpg ~ acceleration, data = auto_mpg)
res_accel <- residuals(fit_accel_ols)

cat("Acceleration residual diagnostics:\n")
cat("  Skewness (gamma3):", round(pmm_skewness(res_accel), 3), "\n")
cat("  Kurtosis (gamma4):", round(pmm_kurtosis(res_accel) - 3, 3), "\n")

# Dispatch
dispatch_accel <- pmm_dispatch(res_accel)

# Compare
fit_accel_pmm2 <- lm_pmm2(mpg ~ acceleration, data = auto_mpg, na.action = na.omit)
vf2_accel <- pmm2_variance_factor(fit_accel_pmm2@m2, fit_accel_pmm2@m3, fit_accel_pmm2@m4)
cat("\nPMM2 g2 =", vf2_accel$g2,
    " (improvement:", round((1 - vf2_accel$g2) * 100, 1), "%)\n")

## ----multiple_regression------------------------------------------------------
n <- 300
x1 <- rnorm(n)
x2 <- runif(n, -2, 2)
x3 <- rnorm(n, 3, 1)

# True parameters
beta <- c(5, 1.5, -0.8, 2.0)

# Beta-symmetric errors (platykurtic, gamma4 ≈ -0.86)
errors_beta <- (rbeta(n, 3, 3) - 0.5) * sqrt(12 * 3)

y_multi <- beta[1] + beta[2]*x1 + beta[3]*x2 + beta[4]*x3 + errors_beta
dat_multi <- data.frame(y = y_multi, x1 = x1, x2 = x2, x3 = x3)

# Fit
fit_multi_ols  <- lm(y ~ x1 + x2 + x3, data = dat_multi)
fit_multi_pmm3 <- lm_pmm3(y ~ x1 + x2 + x3, data = dat_multi)

comparison_multi <- data.frame(
  Parameter = c("Intercept", "x1", "x2", "x3"),
  True  = beta,
  OLS   = coef(fit_multi_ols),
  PMM3  = coef(fit_multi_pmm3),
  Diff  = coef(fit_multi_pmm3) - coef(fit_multi_ols)
)
print(comparison_multi, row.names = FALSE, digits = 4)

cat("\ng3 =", fit_multi_pmm3@g_coefficient,
    " (improvement:", round((1 - fit_multi_pmm3@g_coefficient) * 100, 1), "%)\n")

## ----no_intercept-------------------------------------------------------------
# Model through the origin
y_noint <- 3 * x1 + errors_beta[1:n]
dat_noint <- data.frame(y = y_noint, x1 = x1)

fit_noint <- lm_pmm3(y ~ x1 - 1, data = dat_noint)
cat("True slope: 3\n")
cat("PMM3 estimate:", coef(fit_noint), "\n")
cat("Converged:", fit_noint@convergence, "\n")

## ----interpret_summary--------------------------------------------------------
summary(fit_pmm3)

## ----adaptive_mode------------------------------------------------------------
fit_fixed    <- lm_pmm3(y ~ x, data = dat, adaptive = FALSE)
fit_adaptive <- lm_pmm3(y ~ x, data = dat, adaptive = TRUE)

comparison_adapt <- data.frame(
  Mode       = c("Fixed kappa", "Adaptive kappa"),
  Intercept  = c(coef(fit_fixed)[1], coef(fit_adaptive)[1]),
  Slope      = c(coef(fit_fixed)[2], coef(fit_adaptive)[2]),
  Iterations = c(fit_fixed@iterations, fit_adaptive@iterations)
)
print(comparison_adapt, row.names = FALSE, digits = 5)

## ----diagnostics_full, fig.height=6-------------------------------------------
plot(fit_pmm3)

## ----dispatch_utility---------------------------------------------------------
# Symmetric platykurtic errors → PMM3
pmm_dispatch(runif(500, -1, 1))

# Asymmetric errors → PMM2
pmm_dispatch(rexp(500) - 1)

# Gaussian errors → OLS
pmm_dispatch(rnorm(500))

## ----symmetry_test------------------------------------------------------------
test_symmetry(residuals(fit_ols))

## ----moments_utility----------------------------------------------------------
mom <- compute_moments_pmm3(residuals(fit_ols))
cat("m2 =", mom$m2, "\n")
cat("m4 =", mom$m4, "\n")
cat("m6 =", mom$m6, "\n")
cat("gamma4 =", mom$gamma4, "\n")
cat("kappa =", mom$kappa, "\n")
cat("g3 =", mom$g3, "\n")

## ----gamma6_utility-----------------------------------------------------------
pmm_gamma6(residuals(fit_ols))

## ----variance_factor----------------------------------------------------------
pmm3_variance_factor(mom$m2, mom$m4, mom$m6)

