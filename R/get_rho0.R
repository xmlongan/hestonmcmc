#' Get the sample correlation of the Brownian increments
#'
#' \eqn{\rho_0} is the sample correlation of the Brownian increments of
#' Heston Stochastic Volatility Model
#' \deqn{y_n = \mu h - \frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y}
#' \deqn{v_n - v_{n-1} = k(\theta - v_{n-1})h +\sigma_v\sqrt{v_{n-1}h}\epsilon_n^v}
#' \eqn{\epsilon_n^y} and \eqn{\epsilon_n^v} are first computed from return and
#' volatility samples.
#' \deqn{\epsilon_n^y = \frac{y_n-\mu h + \frac{1}{2}v_{n-1}h}{\sqrt{v_{n-1}h}}
#' ,\qquad
#' \epsilon_n^v=\frac{v_n-v_{n-1}-k(\theta - v_{n-1})h}{\sigma_v\sqrt{v_{n-1}h}}.}
#' Then compute \eqn{\rho_0} as the sample correlation.
#'
#' @param v vector of volatility, \[\eqn{v_0,v_1,\cdots,v_N}\].
#' @param y vector of returns, \[\eqn{y_1,\cdots,y_N}\].
#' @param mu parameter \eqn{\mu}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}, the long-run mean of volatility.
#' @param sigma_v parameter \eqn{\sigma_v}, volatility of volatility.
#' @param h time unit.
#'
#' @return A scale value, \eqn{\rho_0}.
#' @export
#'
#' @examples
#' v = rep(0.25,5)
#' y = rep(0.125,4)
#' mu = 0.125; k = 0.1; theta = 0.25; sigma_v = 0.1
#' h = 1
#' get_rho0(v,y,mu,k,theta,sigma_v,h)
get_rho0 <- function(v,y,mu,k,theta,sigma_v,h) {
  eps_y = rep(0,length(y))
  eps_v = rep(0,length(y))
  for (n in 1:length(y)) {
    eps_y[n] = (y[n] - mu*h + v[n]*h/2)/(sqrt(v[n]*h))
    eps_v[n] = (v[n+1] - v[n] - k*(theta-v[n])*h) / (sigma_v*sqrt(v[n]*h))
  }
  if (stats::sd(eps_y)==0 || stats::sd(eps_v)==0) {
    rho0 = 0
  } else {
    rho0 = stats::cor(eps_y, eps_v)
  }
  if (is.na(rho0)) {
    print(paste('rho0:',rho0))
  }
  return(rho0)
}
