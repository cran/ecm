#'FRED data on Federal Reserve rates and the economy
#'
#'A dataset containing the monthly federal funds rate, unemployment rate, and inflation.
#'@usage data(FedData)
#'@format A data frame with 120 rows and 5 variables:
#'\describe{
#' \item{date}{monthly date}
#' \item{FedFundsRate}{federal funds rate, in percent}
#' \item{UnemploymentRate}{monthly US unemployment rate, in percent}
#' \item{Inflation}{monthly US inflation rate, in percent}
#' \item{GDPgrowth}{quarterly GDP growth, in percent - data have been interpolated and extrapolated to be monthly}
#'}
#'@source \url{https://fred.stlouisfed.org/}
"FedData"