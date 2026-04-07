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
n <- 300
phi_true <- 0.7

# Uniform innovations: symmetric, platykurtic (gamma4 ≈ -1.2)
innov <- runif(n, -sqrt(3), sqrt(3))  # variance = 1

# Generate AR(1) series
ar_data <- arima.sim(n = n, list(ar = phi_true), innov = innov)

# Plot
plot(ar_data, main = "AR(1) Process with Uniform Innovations",
     ylab = "Value", col = "steelblue", lwd = 1.2)

## ----ar1_fit------------------------------------------------------------------
# PMM3 fit
fit_pmm3_ar1 <- ar_pmm3(ar_data, order = 1, include.mean = FALSE)

# MLE fit (classical)
fit_mle_ar1 <- arima(ar_data, order = c(1, 0, 0), include.mean = FALSE)

# Compare
cat("True AR(1) coefficient:", phi_true, "\n")
cat("MLE estimate:          ", coef(fit_mle_ar1)["ar1"], "\n")
cat("PMM3 estimate:         ", coef(fit_pmm3_ar1)["ar1"], "\n")

# PMM3 diagnostics
summary(fit_pmm3_ar1)

## ----ar2_example--------------------------------------------------------------
phi_true <- c(0.5, -0.3)

# Beta-symmetric innovations: symmetric, platykurtic
innov2 <- (rbeta(400, 2, 2) - 0.5) * sqrt(12 * 2)  # scaled for var ≈ 1
ar2_data <- arima.sim(n = 400, list(ar = phi_true), innov = innov2)

fit_pmm3_ar2 <- ar_pmm3(ar2_data, order = 2, include.mean = FALSE)
fit_mle_ar2  <- arima(ar2_data, order = c(2, 0, 0), include.mean = FALSE)

comparison_ar2 <- data.frame(
  Parameter = c("ar1", "ar2"),
  True   = phi_true,
  MLE    = c(coef(fit_mle_ar2)["ar1"], coef(fit_mle_ar2)["ar2"]),
  PMM3   = coef(fit_pmm3_ar2)[c("ar1", "ar2")]
)

print(comparison_ar2, row.names = FALSE)

cat("\nPMM3 g3 =", fit_pmm3_ar2@g_coefficient,
    "(expected variance ratio PMM3/OLS)\n")

## ----ma1_simulation-----------------------------------------------------------
n <- 300
theta_true <- 0.6

# Truncated normal innovations (|x| ≤ 2): symmetric, platykurtic
innov_raw <- rnorm(n * 2)
innov_tn <- innov_raw[abs(innov_raw) <= 2]
innov_tn <- innov_tn[1:n]
innov_tn <- innov_tn / sd(innov_tn)  # standardize

# Generate MA(1) series
ma_data <- arima.sim(n = n, list(ma = theta_true), innov = innov_tn)

# Fit
fit_pmm3_ma1 <- ma_pmm3(ma_data, order = 1, include.mean = FALSE)
fit_mle_ma1  <- arima(ma_data, order = c(0, 0, 1), include.mean = FALSE)

cat("True MA(1) coefficient:", theta_true, "\n")
cat("MLE estimate:          ", coef(fit_mle_ma1)["ma1"], "\n")
cat("PMM3 estimate:         ", coef(fit_pmm3_ma1)["ma1"], "\n")
cat("\nConvergence:", fit_pmm3_ma1@convergence, "\n")
cat("Iterations: ", fit_pmm3_ma1@iterations, "\n")

## ----ma2_example--------------------------------------------------------------
theta_true <- c(0.5, 0.3)

innov_u <- runif(400, -sqrt(3), sqrt(3))
ma2_data <- arima.sim(n = 400, list(ma = theta_true), innov = innov_u)

fit_pmm3_ma2 <- ma_pmm3(ma2_data, order = 2, include.mean = FALSE)
fit_mle_ma2  <- arima(ma2_data, order = c(0, 0, 2), include.mean = FALSE)

comparison_ma2 <- data.frame(
  Parameter = c("ma1", "ma2"),
  True   = theta_true,
  MLE    = c(coef(fit_mle_ma2)["ma1"], coef(fit_mle_ma2)["ma2"]),
  PMM3   = coef(fit_pmm3_ma2)[c("ma1", "ma2")]
)

print(comparison_ma2, row.names = FALSE)

## ----arma11_simulation--------------------------------------------------------
n <- 400
phi_true   <- 0.5
theta_true <- 0.3

# Uniform innovations
innov_arma <- runif(n, -sqrt(3), sqrt(3))

arma_data <- arima.sim(n = n,
                       list(ar = phi_true, ma = theta_true),
                       innov = innov_arma)

# Fit
fit_pmm3_arma <- arma_pmm3(arma_data, order = c(1, 1), include.mean = FALSE)
fit_mle_arma  <- arima(arma_data, order = c(1, 0, 1), include.mean = FALSE)

cat("True: AR =", phi_true, ", MA =", theta_true, "\n\n")

cat("MLE estimates:\n")
print(coef(fit_mle_arma)[c("ar1", "ma1")])

cat("\nPMM3 estimates:\n")
print(coef(fit_pmm3_arma))

summary(fit_pmm3_arma)

## ----arma_diagnostics, fig.height=6-------------------------------------------
plot(fit_pmm3_arma)

## ----arima110_simulation------------------------------------------------------
n <- 300
phi_true <- 0.6

# Generate stationary AR(1) with uniform innovations
innov_arima <- runif(n, -sqrt(3), sqrt(3))
x_stat <- arima.sim(n = n, list(ar = phi_true), innov = innov_arima)

# Integrate (cumulative sum) to create non-stationary series
x_integrated <- cumsum(x_stat)

plot(x_integrated, type = "l",
     main = "ARIMA(1,1,0) Process (Integrated AR(1))",
     ylab = "Value", xlab = "Time", col = "darkgreen", lwd = 1.2)

## ----arima_fit----------------------------------------------------------------
fit_pmm3_arima <- arima_pmm3(x_integrated, order = c(1, 1, 0),
                             include.mean = FALSE)
fit_mle_arima  <- arima(x_integrated, order = c(1, 1, 0),
                        include.mean = FALSE)

cat("True AR(1) coefficient:", phi_true, "\n")
cat("MLE estimate:          ", coef(fit_mle_arima)["ar1"], "\n")
cat("PMM3 estimate:         ", coef(fit_pmm3_arima)["ar1"], "\n")

summary(fit_pmm3_arima)

## ----arima111_example---------------------------------------------------------
phi_true   <- 0.4
theta_true <- 0.3

innov_111 <- runif(350, -sqrt(3), sqrt(3))
x_stat2 <- arima.sim(n = 350, list(ar = phi_true, ma = theta_true),
                     innov = innov_111)
x_int2 <- cumsum(x_stat2)

fit_pmm3_111 <- arima_pmm3(x_int2, order = c(1, 1, 1), include.mean = FALSE)
fit_mle_111  <- arima(x_int2, order = c(1, 1, 1), include.mean = FALSE)

comparison_arima <- data.frame(
  Parameter = c("ar1", "ma1"),
  True   = c(phi_true, theta_true),
  MLE    = c(coef(fit_mle_111)["ar1"], coef(fit_mle_111)["ma1"]),
  PMM3   = coef(fit_pmm3_111)[c("ar1", "ma1")]
)

print(comparison_arima, row.names = FALSE)

## ----forecasting--------------------------------------------------------------
# Forecast 10 steps ahead from the AR(1) model
preds_ar <- predict(fit_pmm3_ar1, n.ahead = 10)

cat("AR(1) PMM3 forecasts (10 steps):\n")
print(round(as.numeric(preds_ar), 4))

# Plot original series with forecasts
n_plot <- 80
plot_data <- tail(as.numeric(ar_data), n_plot)
pred_vals <- as.numeric(preds_ar)

plot(seq_along(plot_data), plot_data, type = "l",
     xlim = c(1, n_plot + 10),
     ylim = range(c(plot_data, pred_vals)),
     main = "AR(1) PMM3: Observed and Forecast",
     xlab = "Time", ylab = "Value",
     col = "steelblue", lwd = 1.5)
lines(n_plot + 1:10, pred_vals, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Observed", "Forecast"),
       col = c("steelblue", "red"), lty = c(1, 2), lwd = c(1.5, 2))

## ----forecast_arma------------------------------------------------------------
# Forecast from the ARMA(1,1) model
preds_arma <- predict(fit_pmm3_arma, n.ahead = 10)

cat("ARMA(1,1) PMM3 forecasts (10 steps):\n")
print(round(as.numeric(preds_arma), 4))

## ----dispatch_example---------------------------------------------------------
# Fit classical model first
fit_classical <- arima(ar_data, order = c(1, 0, 0), include.mean = FALSE)
res_classical <- residuals(fit_classical)

# Dispatch
recommendation <- pmm_dispatch(res_classical)

## ----dispatch_comparison------------------------------------------------------
# Compare all three methods on uniform innovations
cat("Recommended method:", recommendation$method, "\n")
cat("gamma3 =", round(recommendation$gamma3, 3),
    " gamma4 =", round(recommendation$gamma4, 3), "\n")

if (!is.null(recommendation$g3)) {
  cat("g3 =", round(recommendation$g3, 3),
      " (PMM3 variance reduction factor)\n")
}

## ----djia_load----------------------------------------------------------------
data(djia2002)
changes <- na.omit(djia2002$change)
n_obs <- length(changes)

plot(changes, type = "l",
     main = "DJIA Daily Changes (Jul-Dec 2002)",
     ylab = "Change (USD)", xlab = "Trading Day",
     col = "steelblue", lwd = 1.2)
abline(h = 0, col = "red", lty = 2)

## ----djia_analysis------------------------------------------------------------
# Fit a classical AR(1)
fit_djia_mle <- arima(changes, order = c(1, 0, 0))
res_djia <- residuals(fit_djia_mle)

cat("Innovation diagnostics:\n")
cat("  Skewness (gamma3):", round(pmm_skewness(res_djia), 3), "\n")
cat("  Kurtosis (gamma4):", round(pmm_kurtosis(res_djia), 3), "\n")
cat("  Symmetry test:\n")
sym <- test_symmetry(res_djia)
cat("    symmetric:", sym$is_symmetric, "\n")
cat("    message:  ", sym$message, "\n")

# Let pmm_dispatch decide
dispatch <- pmm_dispatch(res_djia)
cat("\nRecommended method:", dispatch$method, "\n")

## ----djia_fit-----------------------------------------------------------------
fit_djia_pmm2 <- ar_pmm2(changes, order = 1)
fit_djia_pmm3 <- ar_pmm3(changes, order = 1)

comparison_djia <- data.frame(
  Method = c("MLE", "PMM2", "PMM3"),
  ar1 = c(coef(fit_djia_mle)["ar1"],
          coef(fit_djia_pmm2)["ar1"],
          coef(fit_djia_pmm3)["ar1"])
)
print(comparison_djia, row.names = FALSE)

## ----oil_load-----------------------------------------------------------------
data(DCOILWTICO)

# Compute log-returns (remove NAs)
prices <- na.omit(DCOILWTICO$DCOILWTICO)
log_returns <- diff(log(prices))
log_returns <- log_returns[is.finite(log_returns)]

plot(log_returns, type = "l",
     main = "WTI Crude Oil: Daily Log-Returns",
     ylab = "Log-return", xlab = "Trading Day",
     col = "darkgreen", lwd = 0.8)
abline(h = 0, col = "red", lty = 2)

## ----oil_analysis-------------------------------------------------------------
# Fit AR(1) to log-returns
fit_oil_mle <- arima(log_returns, order = c(1, 0, 0))
res_oil <- residuals(fit_oil_mle)

cat("Log-return innovation diagnostics:\n")
cat("  Skewness (gamma3):", round(pmm_skewness(res_oil), 3), "\n")
cat("  Kurtosis (gamma4):", round(pmm_kurtosis(res_oil), 3), "\n")

# Dispatch
dispatch_oil <- pmm_dispatch(res_oil)
cat("\nRecommended method:", dispatch_oil$method, "\n")

## ----oil_fit------------------------------------------------------------------
fit_oil_pmm2 <- ar_pmm2(log_returns, order = 1)
fit_oil_pmm3 <- ar_pmm3(log_returns, order = 1)

comparison_oil <- data.frame(
  Method = c("MLE", "PMM2", "PMM3"),
  ar1 = c(coef(fit_oil_mle)["ar1"],
          coef(fit_oil_pmm2)["ar1"],
          coef(fit_oil_pmm3)["ar1"])
)
print(comparison_oil, row.names = FALSE)

vf2_oil <- pmm2_variance_factor(fit_oil_pmm2@m2, fit_oil_pmm2@m3, fit_oil_pmm2@m4)
cat("\nPMM2 g2 =", vf2_oil$g2, "\n")
cat("PMM3 g3 =", fit_oil_pmm3@g_coefficient, "\n")

## ----adaptive_mode------------------------------------------------------------
# Fixed kappa (default)
fit_fixed <- ar_pmm3(ar_data, order = 1, include.mean = FALSE, adaptive = FALSE)

# Adaptive kappa
fit_adapt <- ar_pmm3(ar_data, order = 1, include.mean = FALSE, adaptive = TRUE)

comparison_adaptive <- data.frame(
  Mode   = c("Fixed kappa", "Adaptive kappa"),
  ar1    = c(coef(fit_fixed)["ar1"], coef(fit_adapt)["ar1"]),
  g3     = c(fit_fixed@g_coefficient, fit_adapt@g_coefficient),
  iters  = c(fit_fixed@iterations, fit_adapt@iterations)
)
print(comparison_adaptive, row.names = FALSE)

## ----aic_comparison-----------------------------------------------------------
# Fit several AR orders with PMM3
aic_vals <- sapply(1:5, function(p) {
  fit <- ar_pmm3(ar_data, order = p, include.mean = FALSE)
  AIC(fit)
})

aic_df <- data.frame(Order = 1:5, AIC = round(aic_vals, 2))
print(aic_df, row.names = FALSE)

cat("\nBest AR order by AIC:", which.min(aic_vals), "\n")

