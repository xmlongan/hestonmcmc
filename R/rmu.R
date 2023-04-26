#' Update \eqn{\mu}
#'
#' @description 
#' Parameter \eqn{\mu} has a normal posterior if its prior is also a normal. We
#' generate new \eqn{\mu} from its posterior.
#' * `crmu()` implemented in C++, through package Rcpp.
#' * `rmu()` implemented in R.
#'
#' @details 
#' Employ following version Heston SV model
#' \deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}(\rho\epsilon_n^v
#'       +\sqrt{1-\rho^2}\epsilon_n)}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}\epsilon_n^v}
#' We have
#' \deqn{x_n = \mu + \frac{\sqrt{1-\rho^2}\sqrt{v_{n-1}h}}{h}\epsilon_n}
#' \deqn{x_n|v_{n-1},v_n \sim \mathcal{N}(\mu, \frac{(1-\rho^2)v_{n-1}}{h})}
#' where
#' \deqn{x_n\triangleq \frac{y_n+v_{n-1}h/2-\rho\sqrt{v_{n-1}h}\epsilon_n^v}{h}
#'      ,\quad \sqrt{v_{n-1}h}\epsilon_n^v
#'      = \frac{v_n-v_{n-1}-k(\theta-v_{n-1})h}{\sigma_v}.}
#' The posterior
#' \deqn{\begin{matrix}
#'       P(\mu|v_{0:N},y_{1:N})
#'       &\propto& P(y_{1:N})\cdot P(v_{0:N})\cdot P(\mu)\\
#'       &\propto& \prod_{n=1}^N P(y_n|v_{n-1},v_n) \cdot P(\mu)
#'       \end{matrix}}
#' is a normal when the prior \eqn{P(\mu) = \mathcal{N}(\mu_0,\sigma_0^2)}.
#' The posterior parameters are updated iteratively as
#' \deqn{\sigma_{post}^{2} = \left(\frac{1}{\sigma_0^2} +
#'                             \frac{1}{(1-\rho^2)v_{n-1}/h} \right)^{-1},}
#' \deqn{\mu_{post} = \sigma_{post}^{2}
#'                       \left(\frac{\mu_0}{\sigma_0^2}+
#'                             \frac{x_n}{(1-\rho^2)v_{n-1}/h}\right).}
#' * n=1, compute \eqn{\mu_{post}, \sigma_{post}^2} based on \eqn{v_{n-1},v_n}
#'  and \eqn{\mu_0, \sigma_0^2}.
#' * update prior parameters as
#' \eqn{\mu_0=\mu_{post}, \sigma_0^2=\sigma_{post}^2}.
#' * \eqn{n\ge 2}, repeat above computation until n=N.
#'
#' Generate new \eqn{\mu} using the latest \eqn{\mu_{post}} and
#' \eqn{\sigma_{post}^2} in above repetitve computation.
#'
#' @param prior_mu conjugate prior mean.
#' @param prior_var conjugate prior variance.
#' @param v vector of volatility \[\eqn{v_0, v_1, ..., v_N}\].
#' @param y vector of returns \[\eqn{y_1, ..., y_N}\].
#' @param k parameter k.
#' @param theta parameter \eqn{\theta}.
#' @param sigma_v parameter \eqn{\sigma_v}.
#' @param rho parameter \eqn{\rho}.
#' @param h time unit.
#'
#' @return new value for \eqn{\mu}, a scale value.
#' @export
#'
#' @examples
#' v = rep(0.25,5)
#' y = rep(0.125,4)
#' k = 0.1; theta = 0.25; sigma_v = 0.1; rho = 0
#' h = 1
#' rmu(0,0.1^2,v,y,k,theta,sigma_v,rho,h)
rmu <- function(prior_mu,prior_var,v,y,k,theta,sigma_v,rho,h) {
  N = length(y)
  for (n in 1:N) {
    # v_n <- v_{n-1}
    std_eps_v = ( v[n+1] - v[n] - k*(theta - v[n])*h )/sigma_v # v[1]: v_0
    # y_n <- v_{n-1},v_n
    x_n = ( y[n] + v[n]*h/2 - rho*std_eps_v )/h
    lklhd_var = (1-rho^2)*v[n]/h
    #
    # update parameters, posterior mu and variance
    post_var = 1/(1/prior_var + 1/lklhd_var)
    post_mu = post_var * (prior_mu/prior_var + x_n/lklhd_var)
    #
    prior_mu = post_mu
    prior_var = post_var
  }
  if (post_var < 0) {
    print('post_var < 0 in rmu')
  }
  mu = stats::rnorm(1, mean=post_mu, sd=sqrt(post_var))
  return(mu)
}
