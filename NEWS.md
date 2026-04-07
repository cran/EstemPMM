# EstemPMM News

## Version 0.3.1 (2026-04-06)

### CRAN Resubmission

- Added `CLAUDE.md` and built tarballs to `.Rbuildignore` to eliminate
  non-standard top-level file NOTEs.
- No code changes from 0.3.0; version bump required for CRAN resubmission.

---

## Version 0.3.0 (2026-03-19)

### New Feature: PMM3 for Symmetric Platykurtic Errors

PMM3 (S=3) extends the Polynomial Maximization Method to handle symmetric
error distributions with negative excess kurtosis (platykurtic), such as
uniform, beta-symmetric, and truncated normal errors.

#### New Functions

- **`lm_pmm3()`** - Linear regression estimation using PMM3 (S=3) with
  Newton-Raphson solver. Includes adaptive kappa mode, step-size limiting,
  and divergence guard.
- **`pmm_dispatch()`** - Automatic method selection (OLS / PMM2 / PMM3)
  based on residual cumulant analysis.
- **`compute_moments_pmm3()`** - Compute central moments m2, m4, m6 and
  derived quantities (gamma4, gamma6, g3, kappa) from residuals.
- **`pmm3_variance_factor()`** - Theoretical variance reduction factor
  g3 = 1 - gamma4^2 / (6 + 9*gamma4 + gamma6).
- **`pmm_gamma6()`** - Sixth-order cumulant coefficient.
- **`test_symmetry()`** - Test residual symmetry to guide PMM2 vs PMM3 choice.

#### PMM3 Time Series Functions

- **`ts_pmm3()`** - General PMM3 time series estimation (AR/MA/ARMA/ARIMA).
- **`ar_pmm3()`** - AR model estimation using PMM3.
- **`ma_pmm3()`** - MA model estimation using PMM3.
- **`arma_pmm3()`** - ARMA model estimation using PMM3.
- **`arima_pmm3()`** - ARIMA model estimation using PMM3.

#### New S4 Classes

- **`PMM3fit`** - Standalone class for linear regression (no inheritance
  from BasePMM2) with slots for m2, m4, m6, gamma4, gamma6, g_coefficient,
  and kappa. Full S4 methods: `coef()`, `residuals()`, `fitted()`,
  `predict()`, `summary()`, `plot()`, `AIC()`.
- **`TS3fit`** - Base class for PMM3 time series, with subclasses
  `ARPMM3`, `MAPMM3`, `ARMAPMM3`, `ARIMAPMM3`. Full S4 methods:
  `coef()`, `residuals()`, `fitted()`, `predict()`, `summary()`,
  `plot()`, `AIC()`.

#### Documentation

- New vignette: "PMM3: Linear Regression for Symmetric Platykurtic Errors"
- New vignette: "PMM3 for Time Series: AR, MA, ARMA, and ARIMA Models"
- Updated package-level documentation with PMM3 and method selection sections

---

## Version 0.2.0 (2025-11-20)

### Major Update: Unified PMM2 Architecture

This release represents a significant architectural improvement based on comprehensive research comparing different PMM2 implementation strategies.

#### New Features

- **Unified PMM2 Framework** - Universal PMM2 estimator supporting any nonlinear regression model
  - `pmm2_nonlinear_onestep()` - One-step global correction (default, recommended)
  - `pmm2_nonlinear_iterative()` - Full iterative Newton-Raphson procedure
  - Automatic numerical Jacobian computation via `numDeriv` package
  - Works with AR, MA, ARMA, ARIMA, SAR, SMA, SARIMA models
  
- **Three PMM2 Variants** - New `pmm2_variant` parameter in all time series functions:
  - `"unified_global"` (default) - One-step correction, fast and stable
  - `"unified_iterative"` - Full iterative procedure for maximum accuracy
  - `"linearized"` - Specialized linear approach for MA/SMA models (EstemPMM-style)
  
- **Enhanced Numerical Stability**
  - Optional numerical Jacobian when analytical derivatives unavailable
  - Improved convergence diagnostics
  - Regularization options for ill-conditioned systems

#### Research-Based Improvements

Based on Monte Carlo simulations (R=50, n=200) comparing three approaches:

| Approach | AR(1) | MA(1) | SARIMA | Status |
|----------|-------|-------|---------|--------|
| **Unified Iterative** | -2.9% MSE | -19.9% MSE | **-16.4% MSE** | ✅ **Best overall** |
| **Unified One-step** | -2.2% MSE | **-23.0% MSE** | -15.6% MSE | ✅ **Fastest** |
| **Linearized (MA)** | N/A | -21.6% MSE | N/A | ✅ **MA specialist** |
| Direct Nonlinear | N/A | ❌ **Failed** | ❌ **Failed** | ⛔ **Removed** |

**Key findings:**
- Unified approaches provide consistent 3-23% MSE improvement
- One-step (global) variant offers best speed/accuracy tradeoff
- Linearized approach optimal for pure MA/SMA models

#### Breaking Changes

- **Removed Direct Nonlinear PMM2** - Proved unstable in research (17× worse MSE)
- Previous default behavior preserved with `pmm2_variant = "unified_global"`

#### API Changes

```r
# Old way (still works, uses unified_global by default)
ar_pmm2(y, order = 2)

# New explicit variant selection
ar_pmm2(y, order = 2, pmm2_variant = "unified_iterative")
ma_pmm2(y, order = 1, pmm2_variant = "linearized")  # Best for MA
arima_pmm2(y, order = c(1,0,1), pmm2_variant = "unified_global")  # Default
```

### Documentation

- Updated all function documentation with `pmm2_variant` parameter
- Added comparison table of PMM2 variants to README
- New vignette examples demonstrating variant selection
- Research reports documenting Monte Carlo validation

### Dependencies

- Added `numDeriv` to Suggests for numerical Jacobian computation

### Bug Fixes

- Fixed convergence issues in mixed SARIMA models
- Improved moment estimation for small samples
- Enhanced error messages for degenerate cases

### Performance

- One-step variant: ~50% faster than iterative
- Numerical Jacobian: minimal overhead (<10%) when analytical unavailable
- Memory usage optimized for large time series (n > 1000)

---

## Version 0.1.4 (Development - Superseded by 0.2.0)

### New Features

- **EstemPMM-style PMM2 Estimator for MA/SMA Models** - Advanced parameter estimation for moving average components
  - New `ma_method` argument in `sarima_pmm2()` with options `"mle"` (default) and `"pmm2"`
  - `estpmm_style_ma()` - PMM2 estimator for pure MA(q) models using CSS residuals as fixed regressors
  - `estpmm_style_sma()` - PMM2 estimator for pure SMA(Q) models
  - **`estpmm_style_ma_sma()` - PMM2 estimator for mixed MA+SMA models** ⭐ NEW
  - Full support for MA(q)+SMA(Q) combinations in `sarima_pmm2()` with `ma_method="pmm2"`
  - Expected 20-45% MSE reduction for MA/SMA parameters under asymmetric innovation distributions
  - Implemented in `R/pmm2_ma_estimator.R` module with complete helper functions
  - Comprehensive unit tests (35 total) in `tests/testthat/test-ma-pmm2.R`
  - Addresses limitations identified in Monte Carlo simulations for MA parameter estimation
  - Full backward compatibility - default behavior unchanged

### Bug Fixes

- **Fixed Function Name Conflicts** - Removed obsolete versions of `ma_solve_pmm2`, `ma_compute_innovations`, `sma_compute_innovations`, and `sma_build_design` from `pmm2_ts_main.R` that were overwriting new implementations in `pmm2_ma_estimator.R`
- **Fixed Seasonal Period Validation** - `sarima_pmm2()` now correctly allows `s=1` when no seasonal components (P=0, D=0, Q=0) are specified
- **Corrected ts Object Handling** - MA/SMA estimators now properly convert `ts` objects to numeric vectors before arithmetic operations

## Version 0.1.3 (2025-11-13)

### Documentation

- Expanded both `README.md` and `README_uk.md` with organized function tables, seasonal SAR/SMA workflow examples, and refreshed Monte Carlo efficiency results so new users can discover the seasonal functionality faster.
- Added Part 8 to `vignettes/pmm2_time_series.Rmd`, walking through `sar_pmm2()`/`sma_pmm2()` usage, convergence tips, and practical guidance for seasonal datasets.
- Captured the seasonal-model release summary directly in `NEWS.md`, keeping the changelog aligned with the refreshed documentation.
- Added CRAN-facing housekeeping: refreshed `cran-comments.md`, `CRAN_CHECK_INSTRUCTIONS.md`, and `CRAN_SUBMISSION_CHECKLIST.md`, plus README sections on rebuilding docs, reproducing Monte Carlo studies, and running `R CMD check --as-cran`.


## Version 0.1.2 (2025-11-13)

### New Features

- **Seasonal Autoregressive Models (`sar_pmm2()`)** - Full implementation of SAR(p,P)_s models for seasonal time series
  - Supports arbitrary seasonal periods (e.g., 12 for monthly, 4 for quarterly data)
  - Multiple estimation methods: PMM2, OLS
  - Demonstrated 20-30% variance reduction with asymmetric innovations
  - Full integration with S4 class system (`SARPMM2` class)

    *   Fixed residuals padding in `sar_pmm2`, `sarma_pmm2`, and `sarima_pmm2` to prevent length mismatch errors.
    *   Fixed S4 class definitions to ensure proper method dispatch.
    *   Corrected multiplicative SAR specification in tests.
    *   **Major Improvement**: Enhanced `estpmm_style_ma_sma` to support Multiplicative SARIMA models by including interaction terms in the design matrix. This resolves efficiency issues for mixed MA+SMA models at small sample sizes.
    *   **New Feature**: Extended PMM2 support to full SARIMA models (AR+MA+SAR+SMA) with multiplicative interactions, demonstrating improved efficiency over MLE.


- **Seasonal Moving Average Models (`sma_pmm2()`)** - Complete SMA(Q)_s implementation
  - Flexible seasonal lag specification
  - CSS and PMM2 estimation methods
  - Empirically validated with 500 Monte Carlo replications
  - Achieved 34.1% variance reduction (exceeding theoretical predictions)
  - Robust convergence and computational efficiency

- **Enhanced Comparison Functions**
  - `compare_sar_methods()` - Compare SAR estimation approaches
  - `compare_ts_methods()` - Universal wrapper now supports SAR and SMA models

- **Documentation and Validation**
  - Added comprehensive SAR/SMA documentation in `docs/` directory
  - Monte Carlo validation reports with detailed efficiency metrics
  - Ukrainian language analysis reports
  - Updated both English and Ukrainian READMEs with seasonal models

### Bug Fixes

- **Fixed `predict()` method for PMM2fit class** - The prediction method now correctly handles arbitrary variable names instead of requiring hardcoded "x1", "x2" names. The method now uses general matrix multiplication approach (`X %*% coefficients`) that works with any variable naming convention.
- **Improved coefficient name matching** - Enhanced logic to ensure coefficient names always match design matrix columns, with automatic reordering when necessary.
- **Fixed SAR mean iterations display** - Corrected `sprintf()` call to properly show mean iteration count in comparison output
- **Fixed Seasonal Model Residuals** - `sar_pmm2`, `sarma_pmm2`, and `sarima_pmm2` now correctly pad residuals with zeros (instead of `NA`) to match the original series length, ensuring compatibility with standard diagnostic tools.
- **Fixed S4 Class Definitions** - Reordered class and method definitions in `pmm2_classes.R` to prevent "no definition for class" warnings during package loading.
- **Corrected Multiplicative SAR Specification** - Updated tests to correctly expect 3 coefficients (AR, SAR, Interaction) for multiplicative SAR(1)x(1) models.

### Improvements

- **More robust prediction algorithm** - Simplified prediction code by removing hardcoded special cases and using a unified matrix multiplication approach for all scenarios.
- **Better error messages** - Added clearer error messages when design matrix and coefficients don't match.
- **Enhanced `.gitignore`** - Added `test_results/` directory to version control exclusions

## Version 0.1.1 (2025-10-23)

### Maintenance

- Updated `DESCRIPTION` (latest release date, Suggests list for packages used in the demos).
- Verified the package with `R CMD check --as-cran` (now warning-free after installing `qpdf`).
- Regenerated vignettes (HTML and tangled `.R`) and included them in `inst/doc` for distribution.
- Updated `.Rbuildignore` and `.gitignore`, keeping only files required for CRAN.

## Version 0.1.0 (2025-01-15)

### Initial Release: PMM2 Foundation

**New Features:**
- `lm_pmm2()` - Linear regression estimation using Polynomial Maximization Method (S=2)
- `ar_pmm2()` - Autoregressive (AR) time series modeling with PMM2
- `ma_pmm2()` - Moving Average (MA) time series modeling with PMM2
- `arma_pmm2()` - ARMA time series modeling with PMM2
- `arima_pmm2()` - ARIMA time series modeling with PMM2
- `pmm2_inference()` - Bootstrap inference for linear models
- `ts_pmm2_inference()` - Bootstrap inference for time series models
- Statistical utilities: `pmm_skewness()`, `pmm_kurtosis()`, `compute_moments()`
- Comparison functions: `compare_with_ols()`, `compare_ts_methods()`, `compare_ar_methods()`, `compare_ma_methods()`, `compare_arma_methods()`, `compare_arima_methods()`

**S4 Classes:**
- `PMM2fit` - Results container for linear regression models
- `TS2fit` - Base class for time series results
- `ARPMM2`, `MAPMM2`, `ARMAPMM2`, `ARIMAPMM2` - Specialized time series result classes

**Methods:**
- `summary()` - Model summary statistics
- `coef()` - Extract coefficients
- `fitted()` - Fitted values
- `predict()` - Predictions for new data
- `residuals()` - Model residuals
- `plot()` - Diagnostic plots

**Documentation:**
- Comprehensive Roxygen2 documentation for all exported functions
- README with theoretical background and basic usage examples
- Demonstration script `pmm2_demo_runner.R` showing practical applications

### Package Architecture

**Module Organization:**
- `R/pmm2_main.R` - Primary PMM2 fitting functions
- `R/pmm2_classes.R` - S4 class definitions
- `R/pmm2_utils.R` - Utility functions for moment computation and optimization
- `R/pmm2_ts_design.R` - Time series design matrix construction

**Dependencies:**
- Core: `methods`, `stats`, `graphics`, `utils`
- Optional: `MASS` (for advanced statistical functions, available in Suggests)

**Quality Assurance:**
- Unit tests covering core PMM2 functionality
- Edge case handling for numerical stability
- Convergence diagnostics and warnings

### Known Limitations

- PMM2 only (S=2 order polynomial) - higher orders not yet implemented
- Single-stage estimation (no multi-stage procedures)
- Time series models assume stationarity for AR, MA components
- ARIMA differencing handled via preprocessing, not integrated into core algorithm

### Roadmap

**1.0.0 (Stable API):**
- API stabilization and backward compatibility guarantee
- Seasonal PMM3 models (sar_pmm3, sarima_pmm3)
- Extended performance benchmarks
- Specialized applications (econometrics, biostatistics)

### Citation

If you use EstemPMM in your research, please cite the relevant publications:

**For Linear Regression (lm_pmm2):**
Zabolotnii S., Warsza Z.L., Tkachenko O. (2018) Polynomial Estimation of Linear
Regression Parameters for the Asymmetric PDF of Errors. In: Szewczyk R.,
Zieliński C., Kaliczyńska M. (eds) Automation 2018. AUTOMATION 2018. Advances in
Intelligent Systems and Computing, vol 743. Springer, Cham.
https://doi.org/10.1007/978-3-319-77179-3_75

**For Autoregressive Models (ar_pmm2):**
Zabolotnii S., Tkachenko O., Warsza Z.L. (2022) Application of the Polynomial
Maximization Method for Estimation Parameters of Autoregressive Models with
Asymmetric Innovations. In: Szewczyk R., Zieliński C., Kaliczyńska M. (eds)
Automation 2022. AUTOMATION 2022. Advances in Intelligent Systems and Computing,
vol 1427. Springer, Cham. https://doi.org/10.1007/978-3-031-03502-9_37

**For Moving Average Models (ma_pmm2):**
Zabolotnii S., Tkachenko O., Warsza Z.L. (2023) Polynomial Maximization Method
for Estimation Parameters of Asymmetric Non-gaussian Moving Average Models. In:
Szewczyk R., et al. (eds) Automation 2023. AUTOMATION 2023. Lecture Notes in
Networks and Systems, vol 630. Springer, Cham.

### Technical Notes

**Algorithm Stability:**
- Regularization parameter automatically adjusted for ill-conditioned systems
- Step size limiting prevents divergence in optimization
- Convergence history tracking for diagnostics

**Numerical Considerations:**
- Moment estimation uses robust methods to handle outliers
- Design matrices constructed with numerical stability in mind
- NA/Inf values detected and handled appropriately
