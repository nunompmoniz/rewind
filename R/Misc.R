#' Create Time Delay Embedding
#'
#' Time series forecasting models usually assume the existence of a degree of correlation between successive values of the series. A form of modelling such correlation is to use previous values of the series as predictors of future value(s), in a procedure known as time delay embedding. This process allows the us of standard regression tools on time series forecasting tasks.
#'
#' @param ts the data structure with the time series data
#' @param embed the size of the embed, defaults to 10
#'
#' @return
#' @export
#'
#' @examples
#'
create.data <- function(ts, embed = 10){

  t <- index(ts)[-(1:(embed-1))]
  e <- embed(ts,embed)[, embed:1]
  colnames(e) <- paste('V', 1:embed, sep='')
  d <- xts(e, t)
  as.data.frame(d)

}
