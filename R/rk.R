#' Update \eqn{k}
#'
#' @description 
#' Parameter \eqn{k} has a normal posterior if its prior is also a normal.
#' We generate new \eqn{k} from its posterior.
#' * `crk()` implemented in C++, through package Rcpp.
#' * `rk()` implemented in R.
#'
#' @details 
#' Employ following version Heston SV model
#' \deqn{y_n=\mu h - \frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}
#'      (\rho\epsilon_n^y + \sqrt{1-\rho^2}\epsilon_n)}
#' We have
#' \deqn{x_n = k(\theta-v_{n-1})h +
#'            \sigma_v\sqrt{1-\rho^2}\sqrt{v_{n-1}h}\epsilon_n}
#' \deqn{x_n|v_{n-1},y_n \sim \mathcal{N}(k(\theta-v_{n-1})h,
#'                                         \sigma_v^2(1-\rho^2)v_{n-1}h)}
#' where
#' \deqn{x_n\triangleq v_n-v_{n-1}-\sigma_v\rho\sqrt{v_{n-1}h}\epsilon_n^y
#' , \qquad
#' \sqrt{v_{n-1}h}\epsilon_n^y=y_n-\mu h +\frac{1}{2}v_{n-1}h.}
#' Noting that
#' 1. above representation is robust even for case \eqn{\theta-v_{n-1}=0};
#' 2. `Inf * 0` produces `NaN`, which may happen if we define
#' \eqn{x_n^{'}\triangleq x_n/((\theta-v_{n-1})h)}.
#'
#' The posterior
#' \deqn{\begin{matrix}
#'       P(k|v_{0:N},y_{1:N})
#'       &\propto& P(v_{0:N},y_{1:N}) \cdot P(k)\\
#'       &\propto&\prod_{n=1}^N(P(y_n|v_{n-1})P(v_n|v_{n-1},y_n) )\cdot P(k)\\
#'       &\propto&\prod_{n=1}^NP(v_n|v_{n-1},y_n)\cdot P(k)\\
#'       &\propto&\prod_{n=1}^NP(x_n|v_{n-1},y_n)\cdot P(k)
#'       \end{matrix}}
#' is a normal when the prior \eqn{P(k)=\mathcal{N}(\mu_0,\sigma_0^2)}.
#' The posterior parameters are updated iteratively as
#' \deqn{\sigma_{post}^{2} = \left(\frac{1}{\sigma_0^2} +
#'       \frac{(\theta-v_{n-1})^2h}{\sigma_v^2(1-\rho^2)v_{n-1}} \right)^{-1},}
#' \deqn{\mu_{post} = \sigma_{post}^{2}
#'                       \left(\frac{\mu_0}{\sigma_0^2}+
#'       \frac{x_n(\theta-v_{n-1})}{\sigma_v^2(1-\rho^2)v_{n-1}}\right).}
#' * n=1, compute \eqn{\mu_{post}, \sigma_{post}^2} based on \eqn{v_{n-1},v_n}
#'  and \eqn{\mu_0, \sigma_0^2}.
#' * update prior parameters as
#' \eqn{\mu_0=\mu_{post}, \sigma_0^2=\sigma_{post}^2}.
#' * \eqn{n\ge 2}, repeat above computation until n=N.
#'
#' Generate new \eqn{k} using the latest \eqn{\mu_{post}} and
#' \eqn{\sigma_{post}^2} in above repetitve computation.
#'
#' @param prior_mu conjugate prior mean.
#' @param prior_var conjugate prior variance.
#' @param v vector of volatility \[\eqn{v_0, v_1, ..., v_N}\].
#' @param y vector of returns \[\eqn{y_1, ..., y_N}\].
#' @param mu parameter \eqn{\mu}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma_v parameter \eqn{\sigma_v}.
#' @param rho parameter \eqn{\rho}.
#' @param h time unit.
#'
#' @return new value for \eqn{k}, a scale value.
#' @export
#'
#' @examples
#' v = rep(0.25,5)
#' y = rep(0.125,4)
#' mu = 0; theta = 0.25; sigma_v = 0.1; rho = 0
#' h = 1
#' rk(0.01,1,v,y,mu,theta,sigma_v,rho,h)
rk <- function(prior_mu,prior_var,v,y,mu,theta,sigma_v,rho,h) {
  N = length(y)
  for (n in 1:N) {
    # y_n <- v_{n-1}
    std_eps_y = y[n] - mu*h + v[n]*h/2
    #
    # v_n <- v_{n-1},y_n
    x_n = v[n+1] - v[n] - sigma_v * rho * std_eps_y
    lklhd_var = sigma_v^2 * (1-rho^2) * v[n] * h
    #
    # update parameters, posterior mu and variance
    post_var = 1/(1/prior_var + ((theta-v[n])*h)^2/lklhd_var)
    post_mu = post_var * (prior_mu/prior_var + x_n*(theta-v[n])*h/lklhd_var)
    #
    prior_mu = post_mu
    prior_var = post_var
  }
  k = stats::rnorm(1, mean=post_mu, sd=sqrt(post_var))
  return(k)
}
