#' Simulated returns of Heston Stochastic Volatility Model data
#'
#' A vector of returns with length 1,000, simulated using Euler approximation
#' with true parameter setting:
#' \deqn{(\mu=0.125,k=0.1,\theta=0.25,\sigma_v=0.1,\rho=-0.7), h=1}
#' where \eqn{h} denotes time unit.
#'
#' @format ## `Y_series`
#' A vector with 1,000 elements:
#' \describe{
#'   \item{Y\[n\]}{sample of \eqn{y_n}}
#' }
#' @source Simulated using Euler approximation
"Y_series"
