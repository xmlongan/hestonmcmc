#' Update \eqn{v_i}
#'
#' @description 
#' Generate new intermediate volatility \eqn{v_i} according to its posterior,
#' using Independence Metropolis with a gamma density proposal \eqn{P(\cdot)},
#' i.e., `rgamma(1,shape=2,rate=1/vi_old)`, whose mode equals *vi_old*.
#' Note that the first and last volatility
#' should be generated instead by [rv0()] and [rvN()], respectively.
#' * `crvi()` implemented in C++, through package Rcpp.
#' * `rvi()` implemented in R.
#' 
#' @section Major steps:
#'
#' 1. First, propose a new \eqn{v_i^{(g+1)}} according to gamma distribution,
#' `rgamma(1, shape=2, rate=1/vi_old)`.
#' 2. Then, compute the proposal density ratio
#'    \deqn{\frac{P(v_i^{(g)})}{P(v_i^{(g+1)})}=
#'          \frac{dgamma(v_i^{(g)},shape=2, rate=1/v_i^{(g+1)})}
#'               {dgamma(v_i^{(g+1)},shape=2, rate=1/v_i^{(g)})}.}
#' 3. Last, accept \eqn{v_i^{(g+1)}} with probability
#' \deqn{\alpha(v_i^{(g)},v_i^{(g+1)})
#' =min\left(\frac{\pi(v_i^{(g+1)})}{\pi(v_i^{(g)})}
#'          \times \frac{P(v_i^{(g)})}{P(v_i^{(g+1)})},1\right).}
#'
#' @section Posterior density ratio:
#'
#' Employ following version Heston SV model
#' \deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}(
#'      \rho\epsilon_n^y + \sqrt{1-\rho^2}\epsilon_n)}
#' The posterior
#' \deqn{\begin{matrix}
#'        P(v_i|v_{0:i-1,i+1:N},y_{1:N})
#'        &\propto& P(v_{0:N},y_{1:N})\\
#'        &\propto& P(v_{i+1}|v_i,y_{i+1})\cdot P(y_{i+1}|v_i)
#'        \cdot P(v_i|v_{i-1},y_i)
#'       \end{matrix}}
#' Introduce notation
#' \deqn{\tilde{v}_n\triangleq v_{n-1}+k(\theta-v_{n-1})h+\sigma_v\rho
#'       \sqrt{v_{n-1}h}\epsilon_n^y,\quad
#'       \sqrt{v_{n-1}h}\epsilon_n^y = y_n-\mu h +\frac{1}{2}v_{n-1}h}
#' thus, we have
#' \deqn{v_{i+1}|v_i,y_{i+1}
#'       \sim \mathcal{N}(\tilde{v}_{i+1}, \sigma_v^2(1-\rho^2)v_ih),\quad
#'       P(v_{i+1}|v_i,y_{i+1})
#'       = \frac{1}{\sqrt{2\pi(1-\rho^2)v_ih}\sigma_v}
#'         e^{-\frac{(v_{i+1}-\tilde{v}_{i+1})^2}{2\sigma_v^2(1-\rho^2)v_ih}},}
#' \deqn{y_{i+1}|v_i
#'       \sim \mathcal{N}(\mu h - v_ih/2, v_ih),\quad
#'       P(y_{i+1}|v_i)
#'       = \frac{1}{\sqrt{2\pi v_ih}}
#'         e^{-\frac{(y_{i+1}-\mu h + v_ih/2)^2}{2v_ih}},}
#' \deqn{v_{i}|v_{i-1},y_{i}
#'       \sim \mathcal{N}(\tilde{v}_{i}, \sigma_v^2(1-\rho^2)v_{i-1}h),\quad
#'       P(v_{i}|v_{i-1},y_{i})
#'       = \frac{1}{\sqrt{2\pi(1-\rho^2)v_{i-1}h}\sigma_v}
#'         e^{-\frac{(v_{i}-\tilde{v}_{i})^2}{2\sigma_v^2(1-\rho^2)v_{i-1}h}}.}
#' The posterior density ratio
#' \deqn{\begin{matrix}
#'   \frac{\pi(v_i^{(g+1)})}{\pi(v_i^{(g)})}
#'     &=&\frac{P(v_{i+1}|v_i^{(g+1)},y_{i+1})}{P(v_{i+1}|v_i^{(g)},y_{i+1})}
#'     \times \frac{P(y_{i+1}|v_i^{(g+1)})}{P(y_{i+1}|v_i^{(g)})}
#'     \times \frac{P(v_i^{(g+1)}|v_{i-1},y_i)}{P(v_i^{(g)}|v_{i-1},y_i)}\\
#'   log\left(\frac{\pi(v_i^{(g+1)})}{\pi(v_i^{(g)})}\right)
#'     &=&-log\frac{v_i^{(g+1)}}{v_i^{(g)}}
#'    -\frac{(v_{i+1}-\mu^{(g+1)})^2}{2\sigma_v^2(1-\rho^2)hv_i^{(g+1)}}
#'    +\frac{(v_{i+1}-\mu^{(g)})^2}{2\sigma_v^2(1-\rho^2)hv_i^{(g)}}\\
#'    &&-\frac{(y_{i+1}-\mu h +v_i^{(g+1)}h/2)^2}{2hv_i^{(g+1)}}
#'    +\frac{(y_{i+1}-\mu h +v_i^{(g)}h/2)^2}{2hv_i^{(g)}}\\
#'    &&-\frac{(v_i^{(g+1)}-\tilde{v}_i)^2}{2\sigma_v^2(1-\rho^2)hv_{i-1}}
#'    +\frac{(v_i^{(g)}-\tilde{v}_i)^2}{2\sigma_v^2(1-\rho^2)hv_{i-1}}
#' \end{matrix}}
#' where
#' \deqn{\mu^{(g)}\triangleq v_i^{(g)} + k(\theta-v_i^{(g)})h + \sigma_v\rho
#' \sqrt{v_i^{(g)}h}\epsilon_{i+1}^y}
#' substitute \eqn{g} with \eqn{g+1} to get \eqn{\mu^{(g+1)}}.
#'
#' @param vi_old value of \eqn{v_i} before updating, i.e.,
#' \eqn{v_i^{(g)}}, in the \eqn{g}th iteration.
#' @param vim1 volatility \eqn{v_{i-1}}.
#' @param vip1 volatility \eqn{v_{i+1}}.
#' @param yi return \eqn{y_i}.
#' @param yip1 return \eqn{y_{i+1}}.
#' @param mu parameter \eqn{\mu}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma_v parameter \eqn{\sigma_v}.
#' @param rho parameter \eqn{\rho}.
#' @param h time unit.
#'
#' @return new value of \eqn{v_i}, a scale value.
#' @export
#'
#' @examples
#' v = rep(0.25,5); y = rep(0.125,4)
#' mu = 0.01; k = 0.1; theta = 0.25; sigma_v = 0.1; rho = 0; h = 1
#' i = 1
#' rvi(v[i+1],v[i],v[i+2],y[i],y[i+1],mu,k,theta,sigma_v,rho,h)
rvi <- function(vi_old,vim1,vip1,yi,yip1,mu,k,theta,sigma_v,rho,h) {
  # propose new vi
  vi_new = stats::rgamma(1, shape = 2, rate = 1/vi_old)
  #
  ratio_proposal = exp(
    stats::dgamma(vi_old, shape = 2, rate = 1/vi_new, log = TRUE) -
    stats::dgamma(vi_new, shape = 2, rate = 1/vi_old, log = TRUE))
  #
  # v_i <- v_{i-1},y_i
  mu1 = vim1 + k*(theta-vim1)*h + sigma_v*rho*(yi-mu*h+vim1*h/2)
  var1 = sigma_v^2*(1-rho^2)*vim1*h
  ratio1 = exp(- (vi_new - mu1)^2 / (2*var1)
               + (vi_old - mu1)^2 / (2*var1))
  #
  # y_{i+1} <- v_i
  ratio2 = exp(- (yip1 - mu*h + vi_new*h/2)^2 / (2*vi_new*h)
               + (yip1 - mu*h + vi_old*h/2)^2 / (2*vi_old*h))
  #
  # v_{i+1} <- v_i,y_{i+1}
  mu3_new = vi_new + k*(theta-vi_new)*h + sigma_v*rho*(yip1-mu*h+vi_new*h/2)
  mu3_old = vi_old + k*(theta-vi_old)*h + sigma_v*rho*(yip1-mu*h+vi_old*h/2)
  var_new = sigma_v^2 * (1-rho^2) * vi_new*h
  var_old = sigma_v^2 * (1-rho^2) * vi_old*h
  ratio3 = exp(- (vip1 - mu3_new)^2 / (2*var_new)
               + (vip1 - mu3_old)^2 / (2*var_old))
  #
  ratio_post = (vi_old/vi_new) * ratio3 * ratio2 * ratio1
  #
  ratio = ratio_post * ratio_proposal
  #
  p_accept = min(ratio,1)
  U = stats::runif(1)
  if (U < p_accept) {
    vi = vi_new
  } else {
    vi = vi_old
  }
  return(vi)
}
