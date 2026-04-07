# EstemPMM: Polynomial Maximization Method for Non-Gaussian Regression

[![CRAN status](https://www.r-pkg.org/badges/version/EstemPMM)](https://cran.r-project.org/package=EstemPMM)
[![R](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue)](https://cran.r-project.org/)
[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Overview

**EstemPMM** implements the **Polynomial Maximization Method (PMM)** for parameter
estimation in linear regression and time series models when the error distribution
deviates from normality. PMM achieves lower variance estimates than OLS by exploiting
higher-order moments of the residual distribution.

The package provides two complementary estimators:

| Method | Targets | Key quantity | Typical gain |
|--------|---------|--------------|-------------|
| **PMM2** (S=2) | Asymmetric errors ($\|\gamma_3\| > 0.3$) | $g_2 = 1 - \gamma_3^2/(2+\gamma_4)$ | 10–60% variance reduction |
| **PMM3** (S=3) | Symmetric platykurtic errors ($\gamma_4 < -0.7$) | $g_3 = 1 - \gamma_4^2/(6+9\gamma_4+\gamma_6)$ | 10–40% variance reduction |

Use `pmm_dispatch()` to automatically choose between OLS, PMM2, and PMM3 based on
residual cumulants.

---

## Installation

```r
# Released version from CRAN
install.packages("EstemPMM")

# Development version from GitHub
devtools::install_github("SZabolotnii/EstemPMM")
```

---

## Quick Start

### Automatic Method Selection

```r
library(EstemPMM)

# Fit OLS first, then let pmm_dispatch() choose the best method
fit_ols <- lm(mpg ~ acceleration, data = auto_mpg)
pmm_dispatch(residuals(fit_ols))
# -> Recommends PMM2 (skewed residuals, gamma3 ≈ 0.5)

fit_ols2 <- lm(mpg ~ horsepower + I(horsepower^2), data = na.omit(auto_mpg))
pmm_dispatch(residuals(fit_ols2))
# -> Recommends PMM3 (symmetric platykurtic residuals, gamma4 ≈ -1.3)
```

### PMM2 — Asymmetric Errors

```r
library(EstemPMM)
set.seed(42)

# Simulate regression with skewed (gamma-distributed) errors
n <- 200
x <- rnorm(n)
y <- 2 + 1.5 * x + (rgamma(n, shape = 2, rate = 2) - 1)

# Fit PMM2 and compare with OLS
fit_pmm2 <- lm_pmm2(y ~ x, data = data.frame(y, x))
summary(fit_pmm2)
# g2 < 1 indicates variance reduction over OLS
```

### PMM3 — Symmetric Platykurtic Errors

```r
# Simulate regression with uniform (platykurtic) errors
y3 <- 2 + 1.5 * x + runif(n, -sqrt(3), sqrt(3))

# Fit PMM3 and compare with OLS
fit_pmm3 <- lm_pmm3(y3 ~ x, data = data.frame(y3, x))
summary(fit_pmm3)
# g3 ≈ 0.64 -> 36% variance reduction for uniform errors
```

### Time Series

```r
# AR(1) with asymmetric innovations -> PMM2
y_ar <- arima.sim(n = 300, list(ar = 0.7), innov = rgamma(300, 2, 2) - 1)
fit_ar_pmm2 <- ar_pmm2(y_ar, order = 1)

# AR(1) with uniform innovations -> PMM3
y_ar3 <- arima.sim(n = 300, list(ar = 0.7), innov = runif(300, -sqrt(3), sqrt(3)))
fit_ar_pmm3 <- ar_pmm3(y_ar3, order = 1)

# ARIMA(1,1,1) with PMM3 (uses nonlinear solver)
y_int <- cumsum(arima.sim(n = 300, list(ar = 0.5, ma = 0.3),
                           innov = runif(300, -sqrt(3), sqrt(3))))
fit_arima_pmm3 <- arima_pmm3(y_int, order = c(1, 1, 1))
```

---

## Functions

### Linear Regression

| Function | Description |
|----------|-------------|
| `lm_pmm2()` | PMM2 linear regression (asymmetric errors) |
| `lm_pmm3()` | PMM3 linear regression (symmetric platykurtic errors) |
| `pmm_dispatch()` | Automatic OLS / PMM2 / PMM3 selection |
| `compare_with_ols()` | Side-by-side comparison with OLS |
| `pmm2_inference()` | Bootstrap inference for linear PMM2 models |

### PMM2 Time Series

| Function | Model |
|----------|-------|
| `ar_pmm2()` | AR(p) |
| `ma_pmm2()` | MA(q) |
| `arma_pmm2()` | ARMA(p, q) |
| `arima_pmm2()` | ARIMA(p, d, q) |
| `sar_pmm2()` | Seasonal AR — SAR(p, P)_s |
| `sma_pmm2()` | Seasonal MA — SMA(Q)_s |
| `sarma_pmm2()` | SARMA(p,q)×(P,Q)_s |
| `sarima_pmm2()` | SARIMA(p,d,q)×(P,D,Q)_s |
| `ts_pmm2()` | Universal wrapper |
| `ts_pmm2_inference()` | Bootstrap inference for TS models |

### PMM3 Time Series

| Function | Model |
|----------|-------|
| `ar_pmm3()` | AR(p) |
| `ma_pmm3()` | MA(q) |
| `arma_pmm3()` | ARMA(p, q) |
| `arima_pmm3()` | ARIMA(p, d, q) — nonlinear solver |
| `ts_pmm3()` | Universal wrapper |

### Utilities

| Function | Description |
|----------|-------------|
| `pmm_skewness()` | Skewness coefficient $\gamma_3$ |
| `pmm_kurtosis()` | Excess kurtosis $\gamma_4$ |
| `pmm_gamma6()` | Sixth-order cumulant $\gamma_6$ |
| `test_symmetry()` | Formal symmetry test for residuals |
| `compute_moments()` | Moments m2, m3, m4 |
| `compute_moments_pmm3()` | Moments m2, m4, m6 and PMM3 quantities |
| `pmm2_variance_factor()` | Theoretical $g_2$ factor |
| `pmm3_variance_factor()` | Theoretical $g_3$ factor |

---

## Datasets

| Dataset | Description | N | Use case |
|---------|-------------|---|---------|
| `DCOILWTICO` | WTI crude oil daily prices (FRED) | ~9000 | PMM2 time series |
| `auto_mpg` | UCI Auto MPG vehicle data | 398 | PMM2 and PMM3 linear regression |
| `djia2002` | DJIA daily data, Jul–Dec 2002 | 127 | PMM2 AR(1), published example |

---

## Vignettes

| Vignette | Topic |
|----------|-------|
| `vignette("pmm2-introduction")` | PMM2 linear regression: theory, examples, comparison with OLS |
| `vignette("pmm2-time-series")` | PMM2 for AR/MA/ARMA/ARIMA/SAR/SMA models |
| `vignette("bootstrap-inference")` | Bootstrap confidence intervals for PMM2 |
| `vignette("pmm3-symmetric-errors")` | PMM3 linear regression: uniform/beta/truncated-normal errors |
| `vignette("pmm3-time-series")` | PMM3 for AR/MA/ARMA/ARIMA time series |

---

## S4 Class Hierarchy

```
PMM2fit                      PMM3fit
BasePMM2                     TS3fit
  └─ TS2fit                    ├─ ARPMM3
       ├─ ARPMM2               ├─ MAPMM3
       ├─ MAPMM2               ├─ ARMAPMM3
       ├─ ARMAPMM2             └─ ARIMAPMM3
       ├─ ARIMAPMM2
       ├─ SARPMM2
       ├─ SMAPMM2
       ├─ SARMAPMM2
       └─ SARIMAPMM2
```

All classes provide: `coef()`, `residuals()`, `fitted()`, `predict()`,
`summary()`, `plot()`, `AIC()`.

---

## Efficiency Summary

### PMM2 — Monte Carlo Results (asymmetric innovations)

| Model | Innovation | Variance Reduction |
|-------|------------|-------------------|
| Linear regression | Gamma (shape=2) | 35–57% |
| AR(1) | Exponential | 15–20% |
| MA(1) | Chi-squared (df=3) | 15–23% |
| ARMA(1,1) | Student-t (df=4) | 30–45% |
| SAR(1,1)_12 | Exponential | up to 62% |
| SARMA(1,0,1,1)_12 | Exponential | up to 59% |

### PMM3 — Monte Carlo Results (platykurtic innovations)

| Model | Innovation | $g_3$ | Variance Reduction |
|-------|------------|-------|-------------------|
| Linear regression | Uniform | 0.64 | 36% |
| AR(1) | Uniform | 0.55–0.70 | 30–45% |
| ARIMA(1,1,0) | Uniform | 0.32–0.70 | 30–68% |

---

## PMM2 Variants

Three implementation variants are available for all PMM2 time series functions
via the `pmm2_variant` parameter:

| Variant | Speed | Best for |
|---------|-------|---------|
| `"unified_global"` (default) | Fast | All models, best speed/accuracy tradeoff |
| `"unified_iterative"` | Slower | Maximum accuracy, SARIMA |
| `"linearized"` | Fast | Pure MA/SMA models |

```r
arima_pmm2(y, order = c(1,1,1), pmm2_variant = "unified_iterative")
ma_pmm2(y, order = 1, pmm2_variant = "linearized")
```

---

## References

**Foundational theory:**
Kunchenko, Y.P., Lega, Y.G. (1992). *Estimation of Random Variable Parameters
by the Polynomial Maximization Method*. Naukova Dumka, Kyiv.

**Linear regression (PMM2):**
Zabolotnii S., Warsza Z.L., Tkachenko O. (2018). Polynomial Estimation of Linear
Regression Parameters for the Asymmetric PDF of Errors. *Automation 2018*, AISC
vol. 743. Springer. <https://doi.org/10.1007/978-3-319-77179-3_75>

**Autoregressive models (PMM2):**
Zabolotnii S., Tkachenko O., Warsza Z.L. (2022). Application of the PMM for
Estimation Parameters of Autoregressive Models with Asymmetric Innovations.
*Automation 2022*, AISC vol. 1427. Springer.
<https://doi.org/10.1007/978-3-031-03502-9_37>

**Moving average models (PMM2):**
Zabolotnii S., Tkachenko O., Warsza Z.L. (2023). PMM for Estimation Parameters
of Asymmetric Non-Gaussian Moving Average Models. *Automation 2023*, LNNS
vol. 630. Springer. <https://doi.org/10.1007/978-3-031-25844-2_21>

---

## Author

**Serhii Zabolotnii** — Cherkasy State Business College
ORCID: [0000-0003-0242-2234](https://orcid.org/0000-0003-0242-2234)

Bug reports and feature requests: <https://github.com/SZabolotnii/EstemPMM/issues>

## License

GPL-3 © Serhii Zabolotnii
