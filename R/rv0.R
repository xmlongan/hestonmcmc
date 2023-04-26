#' Update \eqn{v_0}
#'
#' @description 
#' Generate new first volatility \eqn{v_0} according to its posterior, 
#' using Independence
#' Metropolis with a gamma density proposal \eqn{P(\cdot)}, i.e.,
#' `rgamma(1, shape=2, rate=1/v0_old)`, whoes mode equals *v0_old*. 
#' The prior of \eqn{v_0} is assumed to be an uniform over a very large range.
#' (*Depreciated*: Its prior is
#' a gamma distribution with shape \eqn{2k\theta/\sigma_v^2} and rate
#' \eqn{2k/\sigma_v^2}, *requiring* \eqn{k,\theta,\sigma_v} *at least being
#' positive*. Note that the shape and rate are computed and remain
#' unchanged, using the initial values of the parameters.) 
#' * `crv0()` implemented in C++, through package Rcpp.
#' * `rv0()` implemented in R.
#'
#' @section Major steps:
#'
#' 1. First, propose a new \eqn{v_0^{(g+1)}} according to gamma distribution,
#' `rgamma(1, shape=2, rate=1/v0_old)`.
#' 2. Then, compute the proposal density ratio
#'    \deqn{\frac{P(v_0^{(g)})}{P(v_0^{(g+1)})}=
#'          \frac{dgamma(v_0^{(g)},shape=2,rate=1/v_0^{(g+1)})}
#'               {dgamma(v_0^{(g+1)},shape=2,rate=1/v_0^{(g)})}.}
#' 3. Last, accept \eqn{v_0^{(g+1)}} with probability
#' \deqn{\alpha(v_0^{(g)},v_0^{(g+1)})
#' =min\left(\frac{\pi(v_0^{(g+1)})}{\pi(v_0^{(g)})}
#'          \times \frac{P(v_0^{(g)})}{P(v_0^{(g+1)})},1\right).}
#'
#' @section Posterior density ratio:
#'
#' Employ following version Heston SV model
#' \deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}(
#'      \rho\epsilon_n^y + \sqrt{1-\rho^2}\epsilon_n)}
#' The posterior
#' \deqn{\begin{matrix}
#'        P(v_0|v_{1:N},y_{1:N})
#'        &\propto& P(v_0,v_{1:N},y_{1:N})\\
#'        &\propto& P(v_1|v_0,y_1)\cdot P(y_1|v_0) \cdot P_{prior}(v_0)
#'       \end{matrix}}
#' Introduce notation
#' \deqn{\tilde{v}_1\triangleq v_0+k(\theta-v_0)h+\sigma_v\rho
#'       \sqrt{v_0h}\epsilon_1^y,\quad
#'       \sqrt{v_0h}\epsilon_1^y = y_1-\mu h +\frac{1}{2}v_0h}
#' thus, we have
#' \deqn{v_1|v_0,y_1
#'       \sim \mathcal{N}(\tilde{v}_1, \sigma_v^2(1-\rho^2)v_0h),\quad
#'       P(v_1|v_0,y_1)
#'       = \frac{1}{\sqrt{2\pi(1-\rho^2)v_0h}\sigma_v}
#'         e^{-\frac{(v_1-\tilde{v}_1)^2}{2\sigma_v^2(1-\rho^2)v_0h}}.}
#' The posterior density ratio
#' \deqn{\begin{matrix}
#'   \frac{\pi(v_0^{(g+1)})}{\pi(v_0^{(g)})}
#'     &=&\frac{P(v_1|v_0^{(g+1)},y_1)}{P(v_1|v_0^{(g)},y_1)}\times
#'        \frac{P(y_1|v_0^{(g+1)})}{P(y_1|v_0^{(g)})}\times
#'        \frac{P_{prior}(v_0^{(g+1)})}{P_{prior}(v_0^{(g)})}\\
#'   log\left(\frac{\pi(v_0^{(g+1)})}{\pi(v_0^{(g)})}\right)
#'     &=&-log\frac{v_0^{(g+1)}}{v_0^{(g)}}
#'    -\frac{(v_1-\mu^{(g+1)})^2}{2\sigma_v^2(1-\rho^2)hv_0^{(g+1)}}
#'    +\frac{(v_1-\mu^{(g)})^2}{2\sigma_v^2(1-\rho^2)hv_0^{(g)}}\\
#'    &&-\frac{(y_1-\mu h +v_0^{(g+1)}h/2)^2}{2hv_0^{(g+1)}}
#'    +\frac{(y_1-\mu h +v_0^{(g)}h/2)^2}{2hv_0^{(g)}}
#'    +log\frac{P_{prior}(v_0^{(g+1)})}{P_{prior}(v_0^{(g)})}
#' \end{matrix}}
#' where
#' \deqn{\mu^{(g)}\triangleq v_0^{(g)} + k(\theta-v_0^{(g)})h + \sigma_v\rho
#' \sqrt{v_0^{(g)}h}\epsilon_1^y}
#' substitute \eqn{g} with \eqn{g+1} to get \eqn{\mu^{(g+1)}}.
#'
#' @param v0_old value of \eqn{v_0} before updating, i.e.,
#' \eqn{v_0^{(g)}}, in the \eqn{g}th iteration.
#' @param prior_shape parameter shape of \eqn{v_0}'s prior gamma distribution
#' \eqn{P_{prior}(\cdot)}.
#' @param prior_rate parameter rate of \eqn{v_0}'s prior gamma distribution
#' \eqn{P_{prior}(\cdot)}.
#' @param v1 value of \eqn{v_1}.
#' @param y1 value of \eqn{y_1}.
#' @param mu parameter \eqn{\mu}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma_v parameter \eqn{\sigma_v}.
#' @param rho parameter \eqn{\rho}.
#' @param h time unit.
#'
#' @return new \eqn{v_0}, a scale value.
#' @export
#'
#' @examples
#' v1 = 0.2; y1 = 0
#' mu = 0.01; k = 0.1; theta = 0.25; sigma_v = 0.25; rho = 0; h = 1
#' prior_shape = 2*k*theta/sigma_v^2
#' prior_rate = 2*k/sigma_v^2
#' rv0(0.1,prior_shape,prior_rate,v1,y1,mu,k,theta,sigma_v,rho,h)
rv0 <- function(v0_old,prior_shape,prior_rate,v1,y1,mu,k,theta,sigma_v,rho,h) {
  # propose new v0
  v0_new = stats::rgamma(1, shape=2, rate=1/v0_old) # mode = v0_old
  ratio_proposal = exp(
    stats::dgamma(v0_old, shape = 2, rate = 1/v0_new, log = TRUE) -
    stats::dgamma(v0_new, shape = 2, rate = 1/v0_old, log = TRUE))
  #
  # prior density ratio
  # ratio_prior = exp(stats::dgamma(v0_old, prior_shape, prior_rate, log=TRUE) -
  #                   stats::dgamma(v0_new, prior_shape, prior_rate, log=TRUE))
  # above produce error
  ratio_prior = 1 # assuming uniform distribution for the prior
  #
  # y_1 <- v_0
  ratio1 = exp(-(y1 - mu*h + v0_new*h/2)^2 / (2*v0_new*h)
               +(y1 - mu*h + v0_old*h/2)^2 / (2*v0_old*h))
  #
  # v_1 <- v_0,y_1
  m_new = v0_new + k*(theta-v0_new) * h + sigma_v * rho * (y1-mu*h + v0_new*h/2)
  m_old = v0_old + k*(theta-v0_old) * h + sigma_v * rho * (y1-mu*h + v0_old*h/2)
  v_new = sigma_v^2 * (1-rho^2) * v0_new * h
  v_old = sigma_v^2 * (1-rho^2) * v0_old * h
  ratio2 = exp(-(v1 - m_new)^2 / (2*v_new)
               +(v1 - m_old)^2 / (2*v_old))
  #
  ratio_post = (v0_old/v0_new) * ratio2 * ratio1 * ratio_prior
  #
  ratio = ratio_post * ratio_proposal
  if(is.na(ratio)) {
    print(paste0('v0_old: ',v0_old,'; v0_new:',v0_new))
    print(paste0('ratio2: ',ratio2))
    print(paste0('ratio1: ',ratio1))
    print(paste0('ratio_prior: ',ratio_prior))
    print(paste0('ratio_proposal: ',ratio_proposal))
  }
  #
  p_accept = min(ratio,1)
  U = stats::runif(1)
  if (U < p_accept) {
    v0 = v0_new
  } else {
    v0 = v0_old
  }
  return(v0)
}
