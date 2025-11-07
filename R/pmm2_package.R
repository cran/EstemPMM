# pmm2_package.R - Main package file with dependencies and imports

#' @importFrom methods is new slotNames
#' @importFrom graphics abline hist legend lines par plot
#' @importFrom stats acf aggregate arima as.formula cov dnorm lm lm.fit lowess model.frame model.matrix model.response na.fail na.pass pnorm qqline qqnorm quantile arima.sim diffinv rnorm sd
#' @importFrom utils head tail
NULL

#' EstemPMM: Polynomial Maximization Method for Robust Regression and Time Series
#'
#' The EstemPMM package provides robust methods for estimating parameters of linear
#' models and time series models that are robust to non-Gaussian errors.
#'
#' @section Linear Regression Functions:
#'
#' \code{\link{lm_pmm2}} - Fit linear models using PMM2
#'
#' \code{\link{compare_with_ols}} - Compare PMM2 with OLS
#'
#' @section Time Series Functions:
#'
#' \code{\link{ts_pmm2}} - General function for fitting time series models using PMM2
#'
#' \code{\link{ar_pmm2}} - Fit AR models
#'
#' \code{\link{ma_pmm2}} - Fit MA models
#'
#' \code{\link{arma_pmm2}} - Fit ARMA models
#'
#' \code{\link{arima_pmm2}} - Fit ARIMA models
#'
#' \code{\link{compare_ts_methods}} - Compare PMM2 with classical methods
#'
#' @section Statistical Inference:
#'
#' \code{\link{pmm2_inference}} - Bootstrap inference for linear models
#'
#' \code{\link{ts_pmm2_inference}} - Bootstrap inference for time series models
#'
#' @section Utilities:
#'
#' \code{\link{pmm_skewness}} - Compute skewness
#'
#' \code{\link{pmm_kurtosis}} - Compute kurtosis
#'
#' \code{\link{compute_moments}} - Compute moments and cumulants
#' @keywords internal
"_PACKAGE"

#' S4 Class PMM2fit
#'
#' Class for storing results of linear model estimation using PMM2
#'
#' @section Slots:
#' \describe{
#'   \item{coefficients}{Estimated coefficients}
#'   \item{residuals}{Final residuals}
#'   \item{m2}{Second central moment}
#'   \item{m3}{Third central moment}
#'   \item{m4}{Fourth central moment}
#'   \item{convergence}{Convergence status}
#'   \item{iterations}{Number of iterations performed}
#'   \item{call}{Original call}
#' }
#'
#' @docType class
#' @name PMM2fit-class
NULL

#' S4 Class TS2fit
#'
#' Base class for storing results of time series model estimation using PMM2
#'
#' @section Slots:
#' \describe{
#'   \item{coefficients}{Estimated coefficients}
#'   \item{residuals}{Final residuals}
#'   \item{m2}{Second central moment}
#'   \item{m3}{Third central moment}
#'   \item{m4}{Fourth central moment}
#'   \item{convergence}{Convergence status}
#'   \item{iterations}{Number of iterations performed}
#'   \item{call}{Original call}
#'   \item{model_type}{Model type}
#'   \item{intercept}{Intercept}
#'   \item{original_series}{Original time series}
#'   \item{order}{Model orders}
#' }
#'
#' @docType class
#' @name TS2fit-class
NULL
