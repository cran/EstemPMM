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
