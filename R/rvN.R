#' Update \eqn{v_N}
#'
#' @description 
#' Generate new last volatility \eqn{v_N} according to its posterior,
#' using Independence Metropolis with a gamma density proposal \eqn{P(\cdot)},
#' i.e., `rgamma(1, shape=2, rate=1/vN_old)`, whose mode equals *vN_old*.
#' * `crvN()` implemented in C++, through package Rcpp.
#' * `rvN()` implemented in R.
#'
#' @section Major steps:
#'
#' 1. First, propose a new \eqn{v_N^{(g+1)}} according to gamma distribution,
#' `rgamma(1, shape=2, rate=1/vN_old)`.
#' 2. Then, compute the proposal density ratio
#'    \deqn{\frac{P(v_N^{(g)})}{P(v_N^{(g+1)})}=
#'          \frac{dgamma(v_N^{(g)},shape=2,rate=1/v_N^{(g+1)})}
#'               {dgamma(v_N^{(g+1)},shape=2,rate=1/v_N^{(g)})}.}
#' 3. Last, accept \eqn{v_N^{(g+1)}} with probability
#' \deqn{\alpha(v_N^{(g)},v_N^{(g+1)})
#'       =min\left(\frac{\pi(v_N^{(g+1)})}{\pi(v_N^{(g)})}
#'           \times \frac{P(v_N^{(g)})}{P(v_N^{(g+1)})},1\right).}
#'
#' @section Posterior density ratio:
#'
#' Employ following version Heston SV model
#' \deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}(
#'      \rho\epsilon_n^y + \sqrt{1-\rho^2}\epsilon_n)}
#' The posterior
#' \deqn{\begin{matrix}
#'        P(v_N|v_{0:N-1},y_{1:N})
#'        &\propto& P(v_{0:N},y_{1:N})\\
#'        &\propto& P(v_N|v_{N-1},y_N)
#'       \end{matrix}}
#' Introduce notation
#' \deqn{\tilde{v}_N\triangleq v_{N-1}+k(\theta-v_{N-1})h+\sigma_v\rho
#'       \sqrt{v_{N-1}h}\epsilon_N^y,\quad
#'       \sqrt{v_{N-1}h}\epsilon_N^y = y_N-\mu h +\frac{1}{2}v_{N-1}h}
#' thus, we have
#' \deqn{v_{N}|v_{N-1},y_{N}
#'       \sim \mathcal{N}(\tilde{v}_{N}, \sigma_v^2(1-\rho^2)v_{N-1}h),\quad
#'       P(v_{N}|v_{N-1},y_{N})
#'       = \frac{1}{\sqrt{2\pi(1-\rho^2)v_{N-1}h}\sigma_v}
#'         e^{-\frac{(v_{N}-\tilde{v}_{N})^2}{2\sigma_v^2(1-\rho^2)v_{N-1}h}}.}
#' The posterior density ratio
#' \deqn{\begin{matrix}
#'   \frac{\pi(v_N^{(g+1)})}{\pi(v_N^{(g)})}
#'     &=& \frac{P(v_N^{(g+1)}|v_{N-1},y_N)}{P(v_N^{(g)}|v_{N-1},y_N)}\\
#'   log\left(\frac{\pi(v_N^{(g+1)})}{\pi(v_N^{(g)})}\right)
#'     &=& -\frac{(v_N^{(g+1)}-\tilde{v}_N)^2}{2\sigma_v^2(1-\rho^2)hv_{N-1}}
#'         +\frac{(v_N^{(g)}-\tilde{v}_N)^2}{2\sigma_v^2(1-\rho^2)hv_{N-1}}
#' \end{matrix}}
#'
#' @param vN_old value of \eqn{v_N} before updating, i.e.,
#' \eqn{v_N^{(g)}}, in the \eqn{g}th iteration.
#' @param vNm1 volatility \eqn{v_{N-1}}.
#' @param yN return \eqn{y_N}.
#' @param mu parameter \eqn{\mu}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma_v parameter \eqn{\sigma_v}.
#' @param rho parameter \eqn{\rho}.
#' @param h time unit.
#'
#' @return new value of \eqn{v_N}, a scale value.
#' @export
#'
#' @examples
#' v = rep(0.25,5); y = rep(0.125,4)
#' mu = 0.01; k = 0.1; theta = 0.25; sigma_v = 0.1; rho = 0; h = 1
#' N = length(y)
#' rvN(v[N+1],v,y,mu,k,theta,sigma_v,rho,h)
rvN <- function(vN_old,vNm1,yN,mu,k,theta,sigma_v,rho,h) {
  # propose new vN
  vN_new = stats::rgamma(1, shape = 2, rate = 1/vN_old)
  ratio_proposal = exp(
    stats::dgamma(vN_old, shape = 2, rate = 1/vN_new, log = TRUE) -
    stats::dgamma(vN_new, shape = 2, rate = 1/vN_old, log = TRUE))
  #
  # v_N <- v_{N-1},y_N
  muN = vNm1 + k*(theta-vNm1) * h + sigma_v * rho * (yN - mu*h + vNm1*h/2)
  varN = sigma_v^2 * (1-rho^2) * vNm1 * h
  #
  ratio_post = exp(-(vN_new - muN)^2 / (2*varN)
                   +(vN_old - muN)^2 / (2*varN))
  ratio = ratio_post * ratio_proposal
  #
  p_accept = min(ratio,1)
  U = stats::runif(1)
  if (U < p_accept) {
    vN = vN_new
  } else {
    vN = vN_old
  }
  return(vN)
}
