#' Update \eqn{\rho}
#'
#' @description 
#' Update \eqn{\rho}, using Independence Metropolis with proposal density
#' centered at the sample correlation between the Brownian increments.
#' * `crrho()` implemented in C++, through package Rcpp.
#' * `rrho()` implemented in R.
#'
#' @section Major steps:
#'
#' 1. First, compute the density center \eqn{\rho_0}
#' (\eqn{ = corr(\epsilon^y,\epsilon^v)}) of the proposal \eqn{P(\cdot)}, based
#' on samples \[\eqn{v_0, v_1, ..., v_N}\] and \[\eqn{y_1, ..., y_N}\] where
#' \deqn{\epsilon_n^y = \frac{y_n-\mu h+\frac{1}{2}v_{n-1}h}{\sqrt{v_{n-1}h}},
#' \qquad
#' \epsilon_n^v = \frac{v_n-v_{n-1}-k(\theta-v_{n-1})h}{\sigma_v\sqrt{v_{n-1}h}}.}
#'
#' 2. Then, propose a new \eqn{\rho^{(g+1)}} according to \eqn{P(\cdot)}, see
#' [rrho_proposal()], [drho_proposal()] for details.
#'
#' 3. Last, accept \eqn{\rho^{(g+1)}} with probability
#' \deqn{\alpha(\rho^{(g)},\rho^{(g+1)})
#' =min\left(\frac{\pi(\rho^{(g+1)})}{\pi(\rho^{(g)})}
#'          \times \frac{P(\rho^{(g)})}{P(\rho^{(g+1)})}, 1\right).}
#'
#' @section Posterior density ratio:
#'
#' Employ following version Heston SV model
#' \deqn{y_n=\mu h-\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}(\rho\epsilon_n^v
#'           +\sqrt{1-\rho^2}\epsilon_n)}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}\epsilon_n^v}
#' The posterior of \eqn{\rho}
#' \deqn{\begin{matrix}
#'       P(\rho|v_{0:N},y_{1:N})
#'       &\propto& P(y_{1:N}|v_{0:N})\cdot P(v_{0:N})\cdot P(\rho)\\
#'       &\propto& \prod_{n=1}^{N}P(y_n|v_{n-1},v_n)\cdot P(\rho)\\
#'       &\propto& \prod_{n=1}^{N}P(y_n|v_{n-1},v_n)\\
#'       &\propto& \prod_{n=1}^{N}P(x_n|v_{n-1},v_n)
#' \end{matrix}}
#' where \eqn{x_n\triangleq y_n-\mu h+\frac{1}{2}v_{n-1}h}, \eqn{P(\rho)} denote
#' the prior density, not the proposal density.
#' Note that \eqn{P(\rho)} is uniform(-1,1), thus constant.
#' \deqn{x_n|v_{n-1},v_n \sim \mathcal{N}(\rho\sqrt{v_{n-1}h}\epsilon_n^v,
#'       (1-\rho^2)v_{n-1}h)}
#' where
#' \deqn{\sqrt{v_{n-1}h}\epsilon_n^v =
#'       \frac{v_n-v_{n-1}-k(\theta-v_{n-1})h}{\sigma_v}.}
#' \deqn{P(x_n|v_{n-1},v_n)=\frac{1}{\sqrt{2\pi(1-\rho^2)v_{n-1}h}}
#'       e^{-\frac{(x_n-\rho\sqrt{v_{n-1}h}\epsilon_n^v)^2}{2(1-\rho^2)v_{n-1}h}}.}
#' Then
#' \deqn{\frac{\pi(\rho^{(g+1)})}{\pi(\rho^{(g)})} =
#'       \frac{P(\rho^{(g+1)}|v_{0:N},y_{1:N})}{P(\rho^{(g)}|v_{0:N},y_{1:N})}.}
#' \deqn{log\left(\frac{\pi(\rho^{(g+1)})}{\pi(\rho^{(g)})}\right)
#'       = - \frac{N}{2}log\frac{1-\rho^{(g+1)2}}{1-\rho^{(g)2}}
#'       - \sum_{n=1}^{N}\left[\frac{(x_n-\rho^{(g+1)}\sqrt{v_{n-1}h}\epsilon_n^v)^2}
#'       {2(1-\rho^{(g+1)2})v_{n-1}h}
#'       - \frac{(x_n-\rho^{(g)}\sqrt{v_{n-1}h}\epsilon_n^v)^2}
#'       {2(1-\rho^{(g)2})v_{n-1}h}\right]
#'       .}
#' Note that \eqn{\rho} can't be -1 or 1, otherwise there will produce `Inf`
#' and may produce `NaN`.
#'
#' @param rho_old value of \eqn{\rho} before updating, \eqn{\rho^{(g)}}, in
#' \eqn{g}th iteration.
#' @param v vector of volatility \[\eqn{v_0, v_1, ..., v_N}\].
#' @param y vector of returns \[\eqn{y_1, ..., y_N}\].
#' @param mu parameter \eqn{\mu}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma_v parameter \eqn{\sigma_v}.
#' @param h time unit.
#'
#' @return new value for \eqn{\rho}, a scale value.
#' @seealso [rrho_proposal()], [drho_proposal()] for information about the
#' proposal distribution \eqn{P(\cdot)}.
#' @export
#'
#' @examples
#' rho_old = 0
#' v = rep(0.25,5)
#' y = rep(0.125,4)
#' mu = 0.125; k = 0.1; theta = 0.25; sigma_v = 0.1
#' h = 1
#' rrho(rho_old,v,y,mu,k,theta,sigma_v,h)
rrho <- function(rho_old,v,y,mu,k,theta,sigma_v,h) {
  # compute the sample correlation between the Brownian increments
  rho0 = get_rho0(v,y,mu,k,theta,sigma_v,h)
  # propose a new rho
  rho_new = rrho_proposal(rho0)
  ratio_proposal = drho_proposal(rho_old, rho0)/
                   drho_proposal(rho_new, rho0)
  #
  sum_old = 0; sum_new =0; N = length(y)
  for (n in 1:N) {
    # y_n <- v_{n-1},v_n
    x_n = y[n] - mu*h + v[n]*h/2
    # v_n <- v_{n-1}
    std_eps_v = (v[n+1] - v[n] - k*(theta - v[n])*h ) / sigma_v
    sum_old = sum_old + (x_n - rho_old*std_eps_v)^2 / (2*(1-rho_old^2)*v[n]*h)
    sum_new = sum_new + (x_n - rho_new*std_eps_v)^2 / (2*(1-rho_new^2)*v[n]*h)
  }
  ratio_post = exp( -log((1-rho_new^2)/(1-rho_old^2))*N/2 - sum_new + sum_old )
  #
  ratio = ratio_post * ratio_proposal
  #
  p_accept = min(ratio,1)
  U = stats::runif(1)
  if (U < p_accept) {
    rho = rho_new
  } else {
    rho = rho_old
  }
  return(rho)
}
