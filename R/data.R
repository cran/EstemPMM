#' WTI Crude Oil Prices
#'
#' Daily spot prices for West Texas Intermediate (WTI) crude oil
#' in U.S. dollars per barrel.
#'
#' @format A data frame with observations for each trading day:
#' \describe{
#'   \item{observation_date}{Date of observation in YYYY-MM-DD format}
#'   \item{DCOILWTICO}{Crude Oil Price: West Texas Intermediate (WTI) in USD per barrel}
#' }
#'
#' @source Federal Reserve Economic Data (FRED), Federal Reserve Bank of St. Louis
#'   \url{https://fred.stlouisfed.org/series/DCOILWTICO}
#'
#' @examples
#' data(DCOILWTICO)
#' head(DCOILWTICO)
#' summary(DCOILWTICO$DCOILWTICO)
"DCOILWTICO"


#' Auto MPG Dataset
#'
#' Fuel consumption and vehicle characteristics for 398 automobiles
#' from the 1970s and 1980s. This dataset is used in published PMM research
#' to demonstrate both PMM2 (asymmetric residuals: MPG vs Weight) and
#' PMM3 (symmetric platykurtic residuals: MPG vs Horsepower).
#'
#' @format A data frame with 398 rows and 9 variables:
#' \describe{
#'   \item{mpg}{Miles per gallon (fuel efficiency)}
#'   \item{cylinders}{Number of cylinders (4, 6, or 8)}
#'   \item{displacement}{Engine displacement (cubic inches)}
#'   \item{horsepower}{Engine horsepower (6 missing values)}
#'   \item{weight}{Vehicle weight (pounds)}
#'   \item{acceleration}{Time to accelerate from 0 to 60 mph (seconds)}
#'   \item{model_year}{Model year (70-82, i.e., 1970-1982)}
#'   \item{origin}{Origin (1 = American, 2 = European, 3 = Japanese)}
#'   \item{car_name}{Car model name}
#' }
#'
#' @details
#' Three regression examples from published PMM papers:
#' \itemize{
#'   \item \strong{MPG vs Acceleration} (PMM2, linear): residuals have
#'     gamma3 = 0.49, g2 = 0.86 (Zabolotnii et al., 2018)
#'   \item \strong{MPG vs Weight} (PMM2, quadratic): residuals have
#'     gamma3 = 0.8, g2 = 0.83 (Zabolotnii et al., 2025)
#'   \item \strong{MPG vs Horsepower} (PMM3, quadratic): residuals have
#'     gamma3 ~ 0.2, gamma4 = 1.3, g3 = 0.89 (Zabolotnii et al., 2025)
#' }
#'
#' @source UCI Machine Learning Repository
#'   \url{https://archive.ics.uci.edu/dataset/9/auto+mpg}
#'
#' @references
#' Quinlan, J.R. (1993). Combining Instance-Based and Model-Based Learning.
#' In Proceedings on the Tenth International Conference of Machine Learning,
#' 236-243.
#'
#' Zabolotnii S., Warsza Z.L., Tkachenko O. (2018) Polynomial Estimation
#' of Linear Regression Parameters for the Asymmetric PDF of Errors.
#' Springer AISC, vol 743. \doi{10.1007/978-3-319-77179-3_75}
#'
#' @examples
#' data(auto_mpg)
#' # PMM2 example: MPG vs Acceleration (asymmetric residuals)
#' fit_ols <- lm(mpg ~ acceleration, data = auto_mpg)
#' pmm_skewness(residuals(fit_ols))  # gamma3 ~ 0.5 -> PMM2
#' pmm_dispatch(residuals(fit_ols))
#' fit_pmm2 <- lm_pmm2(mpg ~ acceleration, data = auto_mpg, na.action = na.omit)
#' coef(fit_pmm2)  # compare with coef(fit_ols)
"auto_mpg"


#' Dow Jones Industrial Average Daily Data (July-December 2002)
#'
#' Daily closing prices and changes of the Dow Jones Industrial Average
#' for the second half of 2002. Used in published PMM2 research to
#' demonstrate AR(1) estimation with asymmetric innovations.
#'
#' @format A data frame with 127 rows and 3 variables:
#' \describe{
#'   \item{date}{Trading date}
#'   \item{close}{DJIA closing price (USD)}
#'   \item{change}{Daily change in closing price (first difference;
#'     NA for the first observation)}
#' }
#'
#' @details
#' The daily changes exhibit positive skewness, making this dataset
#' suitable for PMM2 estimation. In the original paper, an AR(1) model
#' fitted to the change series yielded PMM2 coefficient a1 = -0.43
#' versus OLS a1 = -0.49, with g2 = 0.77 (23\% variance reduction).
#'
#' @source Yahoo Finance via the quantmod R package.
#'
#' @references
#' Zabolotnii S., Tkachenko O., Warsza Z.L. (2022) Application of the
#' Polynomial Maximization Method for Estimation Parameters of Autoregressive
#' Models with Asymmetric Innovations. Springer AISC, vol 1427.
#' \doi{10.1007/978-3-031-03502-9_37}
#'
#' @examples
#' data(djia2002)
#' # AR(1) with PMM2
#' changes <- na.omit(djia2002$change)
#' pmm_skewness(changes)  # positive skewness -> PMM2
#' fit <- ar_pmm2(changes, order = 1)
#' summary(fit)
"djia2002"
