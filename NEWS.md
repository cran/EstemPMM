# EstemPMM News

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

### Roadmap for Future Versions

**0.2.0 (PMM3 Ready Architecture):**
- PMM3 implementation (S=3 polynomial methods)
- Refactored base classes supporting method extensibility
- Vignette documentation with practical use cases
- Enhanced bootstrap procedures for small samples
- GitHub Actions CI/CD integration

**0.3.0 (Advanced Methods):**
- Adaptive PMM order selection
- Robust variance estimation
- Model selection criteria (AIC/BIC for PMM)
- Generalized Linear Models (GLM) with PMM

**1.0.0 (Stable API):**
- API stabilization and backward compatibility guarantee
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
