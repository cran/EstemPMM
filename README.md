# EstemPMM: Polynomial Maximization Method for Regression Analysis

[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://cran.r-project.org/)
[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Overview

The EstemPMM package implements the Polynomial Maximization Method (PMM) for estimating linear regression parameters in cases where the error distribution differs from normal, particularly when it has an asymmetric character.

PMM allows obtaining parameter estimates with lower variance compared to the classical Ordinary Least Squares (OLS) method, especially when the error distribution has significant asymmetry.

## Theoretical Background

PMM uses polynomials of degree S for parameter estimation. When S=1, PMM estimates coincide with OLS estimates. When S=2 and there is asymmetry in the error distribution, PMM can provide reduced variance of estimates.

The PMM framework was originally developed by Yu. P. Kunchenko, who formulated the polynomial maximization approach for estimating distribution parameters [Kunchenko & Lega, 1992].

The theoretical coefficient of variance reduction for S=2 is calculated by the formula:

```
g = 1 - c3^2 / (2 + c4)
```

where `c3` is the skewness coefficient and `c4` is the kurtosis coefficient.


## Installation

```r
# Install from GitHub (requires 'devtools' package)
devtools::install_github("SZabolotnii/EstemPMM")
```

## Basic Usage

```r
library(EstemPMM)

# Create data with asymmetric errors
n <- 100
x <- rnorm(n)
errors <- rgamma(n, shape = 2, scale = 1) - 2  # Shift for zero mean
y <- 2 + 1.5 * x + errors
data <- data.frame(x = x, y = y)

# Fit the model using PMM2
fit <- lm_pmm2(y ~ x, data = data)

# Review results
summary(fit)

# Compare with OLS
ols_fit <- lm(y ~ x, data = data)
compare_with_ols(y ~ x, data)
```

## Main Functions

- `lm_pmm2()`: Fit linear regression using PMM for S=2
- `summary()`: Display fitting results
- `predict()`: Make predictions based on a PMM model
- `pmm2_inference()`: Statistical inference through bootstrap
- `compare_with_ols()`: Compare with OLS estimates
- `plot()`: Diagnostic plots for PMM models

## Project Structure

The package consists of several key R files:
- `pmm2_classes.R`: Defines S4 classes for PMM2 fit results
- `pmm2_main.R`: Contains linear PMM2 fitting functions
- `pmm2_ts_main.R`: Provides PMM2 fitting wrappers for time series models
- `pmm2_utils.R`: Provides optimization utilities for PMM2
- `pmm2_common.R`: Hosts shared numerical routines used by PMM2 algorithms
- `pmm2_inference.R`: Implements bootstrap inference for PMM2 fits
- `pmm2_ts_design.R`: Builds design matrices and helpers for time series estimation
- `pmm2_simulation.R`: Contains code for Monte Carlo simulations
- `pmm2_real_data.R`: Applies PMM2 to real-world data (Auto MPG dataset)

## Variance Diagnostics

You can inspect the theoretical skewness, kurtosis, and expected variance reduction
obtained from the PMM2 algorithm using the helper utilities:

```r
library(EstemPMM)

set.seed(42)
x <- rnorm(200)
errors <- rgamma(200, shape = 2, scale = 1) - 2
y <- 1 + 0.7 * x + errors
dat <- data.frame(y, x)

fit <- lm_pmm2(y ~ x, data = dat)

# Retrieve theoretical cumulants and variance ratio g
vf <- pmm2_variance_factor(fit@m2, fit@m3, fit@m4)
vf$g  # Expected Var(PMM2) / Var(OLS)

# Compare variance matrices for OLS and PMM2
vm <- pmm2_variance_matrices(attr(fit, "model_matrix"),
                             fit@m2, fit@m3, fit@m4)
vm$pmm2
```

## Demo Script

The package includes a detailed demonstration script `pmm2_demo_runner.R` that shows:

1. Comparison of PMM2 and OLS on data with different error distributions
2. Monte Carlo simulations for efficiency evaluation
3. Application to real data (Auto MPG dataset)
4. Bootstrap analysis for uncertainty estimation

To run the demonstration, execute:

```r
source("pmm2_demo_runner.R")
all_results <- run_all_simulations()  # For Monte Carlo simulations
results <- apply_to_mpg_data()  # For real data analysis
```

## Results and Efficiency

PMM2 is particularly effective for distributions with high asymmetry:

| Distribution    | Skewness | Kurtosis | Theoretical Improvement | Actual Improvement |
|-----------------|----------|----------|------------------------|-------------------|
| Gamma (a=0.5)   | 2.83     | 12       | 57%                    | ~50%              |
| Exponential     | 2.00     | 6        | 50%                    | ~45%              |
| Gamma (a=2)     | 1.41     | 3        | 40%                    | ~35%              |
| Lognormal       | 1.00     | 1.5      | 29%                    | ~25%              |

## Adaptive Estimation

The package implements an adaptive procedure for PMM estimation:
1. Find initial OLS estimates and calculate residuals
2. Estimate moments and cumulants of the OLS residuals
3. Calculate refined PMM estimates using these moment estimates

This approach doesn't require prior knowledge of the error distribution properties.

## Applications

The method is particularly useful in:
- Economic and financial modeling with asymmetric error distributions
- Biological systems analysis
- Technical measurements with non-Gaussian noise
- Any regression problem where error distributions exhibit significant skewness

## Authors

- Serhii Zabolotnii - Cherkasy State Business College


## Scientific Publications

The Polynomial Maximization Method and its applications are described in the following peer-reviewed publications:

### Foundational Reference
Kunchenko, Y. P., & Lega, Y. G. (1992). *Estimation of Random Variable Parameters by the Polynomial Maximization Method*. Kyiv: Naukova Dumka. 180 pp.

### Linear Regression (PMM2 for lm_pmm2)
Zabolotnii S., Warsza Z.L., Tkachenko O. (2018) Polynomial Estimation of Linear Regression Parameters for the Asymmetric PDF of Errors. In: Szewczyk R., Zieliński C., Kaliczyńska M. (eds) Automation 2018. AUTOMATION 2018. Advances in Intelligent Systems and Computing, vol 743. Springer, Cham. https://doi.org/10.1007/978-3-319-77179-3_75

### Autoregressive Models (PMM2 for ar_pmm2)
Zabolotnii S., Tkachenko O., Warsza Z.L. (2022) Application of the Polynomial Maximization Method for Estimation Parameters of Autoregressive Models with Asymmetric Innovations. In: Szewczyk R., Zieliński C., Kaliczyńska M. (eds) Automation 2022. AUTOMATION 2022. Advances in Intelligent Systems and Computing, vol 1427. Springer, Cham. https://doi.org/10.1007/978-3-031-03502-9_37

### Moving Average Models (PMM2 for ma_pmm2)
Zabolotnii S., Tkachenko O., Warsza Z.L. (2023) Polynomial Maximization Method for Estimation Parameters of Asymmetric Non-gaussian Moving Average Models. In: Szewczyk R., et al. (eds) Automation 2023. AUTOMATION 2023. Lecture Notes in Networks and Systems, vol 630. Springer, Cham. https://doi.org/10.1007/978-3-031-25844-2_21

## License

This project is distributed under the GPL-3 License.
