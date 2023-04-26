#' Update \eqn{\sigma_v}
#'
#' @description 
#' Update \eqn{\sigma_v}, using Independence Metropolis with a gamma
#' distribution as the
#' proposal distribution. Meanwhile, the prior distribution of \eqn{\sigma_v}
#' is also (assumed) a gamma distribution with *prior_shape* and *prior_rate*.
#' Note that these two prior parameters should be chosen carefully to avoid
#' `Inf*0` in the middle of the computation. The prior density shouldn't be
#' close to 0 for the effective range of \eqn{\sigma_v}.
#' * `crsigma_v()` implemented in C++, through package Rcpp.
#' * `rsigma_v()` implemented in R.
#'
#' @section Major steps:
#'
#' 1. First, compute the rate parameter of the proposal distribution
#' \eqn{P(\cdot)} (gamma distribution, the shape parameter defaults to 2). Since
#' \deqn{var(v_n-v_{n-1})\approx\sigma_v^2E[v_{n-1}]h,}
#' we have
#' \deqn{\sigma_v^2\approx \frac{S^2_{v_n-v_{n-1}}}{\bar{v}_{n-1}h}}
#' where \eqn{S^2_{v_n-v_{n-1}}}, \eqn{\bar{v}_{n-1}} are the
#' corresponding sample variance and mean, respectively.
#' If we want the proposal has a mode
#' equals to square root of above approximation, we can set the proposal
#' rate as
#' \deqn{\beta_{proposal} = 1/
#'       \sqrt{S^2_{v_n-v_{n-1}}/(\bar{v}_{n-1}h)}.}
#' Note that mode of gamma distribution equals \eqn{(\alpha -1)/\beta} for
#' shape \eqn{\alpha \ge 1} and rate \eqn{\beta}.
#' **Warning**: \eqn{S_{v_n-v_{n-1}}^2=0} happens in extreme case!
#'
#' 2. Then, propose a new \eqn{\sigma_v^{(g+1)}} according to gamma
#' distribution, `rgamma(1,shape=2,rate)`, where rate = \eqn{\beta_{proposal}}.
#'  The proposal density ratio
#' \deqn{\frac{P(\sigma_v^{(g)})}{P(\sigma_v^{(g+1)})}
#'    =min\left( \frac{dgamma(\sigma_v^{(g)}, shape=2, rate=1/\beta_{proposal})}
#'    {dgamma(\sigma_v^{(g+1)}, shape=2, rate=1/\beta_{proposal})}
#'    \right).}
#'
#' 3. Last, accept \eqn{\sigma_v^{(g+1)}} with probability
#' \deqn{\alpha(\sigma_v^{(g)},\sigma_v^{(g+1)})
#' =min\left(\frac{\pi(\sigma_v^{(g+1)})}{\pi(\sigma_v^{(g)})}
#'          \times \frac{P(\sigma_v^{(g)})}{P(\sigma_v^{(g+1)})}, 1\right).}
#'
#' @section Posterior density ratio:
#'
#' Employ following version Heston SV model
#' \deqn{y_n=\mu h - \frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}
#'       (\rho\epsilon_n^y + \sqrt{1-\rho^2}\epsilon_n)}
#' The posterior of \eqn{\sigma_v}
#' \deqn{\begin{matrix}
#'   P(\sigma_v|v_{0:N},y_{1:N})
#'   &\propto& P(y_{1:N},v_{0:N})\cdot P_{prior}(\sigma_v)\\
#'   &\propto&\prod_{n=1}^{N}(P(v_n|v_{n-1},y_n)P(y_n|v_{n-1}))
#'    \cdot P_{prior}(\sigma_v)\\
#'   &\propto&\prod_{n=1}^{N}P(v_n|v_{n-1},y_n)\cdot P_{prior}(\sigma_v)\\
#'   &\propto&\prod_{n=1}^{N}P(x_n|v_{n-1},y_n)\cdot P_{prior}(\sigma_v)
#' \end{matrix}}
#' where \eqn{x_n\triangleq v_n-v_{n-1}-k(\theta-v_{n-1})h},
#' \eqn{P_{prior}(\sigma_v)}
#' denotes the prior density of \eqn{\sigma_v}. We have
#' \deqn{x_n|v_{n-1},y_n\sim \mathcal{N}(
#'      \sigma_v\rho\sqrt{v_{n-1}h}\epsilon_n^y,
#'      \sigma_v^2(1-\rho^2)v_{n-1}h)}
#' where
#' \deqn{\sqrt{v_{n-1}h}\epsilon_n^y=y_n-\mu h+\frac{1}{2}v_{n-1}h.}
#' \deqn{P(x_n|v_{n-1},y_n) = \frac{1}{\sigma_v\sqrt{2\pi(1-\rho^2)v_{n-1}h}}
#'       e^{-\frac{(x_n-\sigma_v\rho\sqrt{v_{n-1}h}\epsilon_n^y)^2}
#'       {2\sigma_v^2(1-\rho^2)v_{n-1}h}}.}
#'
#' Then the posterior density ratio
#' \deqn{\frac{\pi(\sigma_v^{(g+1)})}
#'            {\pi(\sigma_v^{(g)})} =
#'       \frac{P(\sigma_v^{(g+1)}|v_{0:N},y_{1:N})}
#'            {P(\sigma_v^{(g)}|v_{0:N},y_{1:N})}.}
#' \deqn{log\left(\frac{\pi(\sigma_v^{(g+1)})}{\pi(\sigma_v^{(g)})}\right)
#'   = - N\cdot log\frac{\sigma_v^{(g+1)}}{\sigma_v^{(g)}}
#'   - \sum_{n=1}^N\left[
#'   \frac{(x_n-\sigma_v^{(g+1)}\rho\sqrt{v_{n-1}h}\epsilon_n^y)^2}
#'   {2\sigma_v^{(g+1)2}(1-\rho^2)v_{n-1}h}
#'   - \frac{(x_n-\sigma_v^{(g)}\rho\sqrt{v_{n-1}h}\epsilon_n^y)^2}
#'   {2\sigma_v^{(g)2}(1-\rho^2)v_{n-1}h}\right]
#'   + log\frac{P_{prior}(\sigma_v^{(g+1)})}{P_{prior}(\sigma_v^{(g)})}.}
#'
#' @param sigma_v_old value of \eqn{\sigma_v} before updating, i.e.,
#' \eqn{\sigma_v^{(g)}}, in \eqn{g}th iteration.
#' @param prior_shape parameter shape of the prior of \eqn{\sigma_v}, (assumed)
#' gamma distribution.
#' @param prior_rate parameter rate of the prior of \eqn{\sigma_v}, (assumed)
#' gamma distribution.
#' @param v vector of volatility \[\eqn{v_0, v_1, ..., v_N}\].
#' @param y vector of returns \[\eqn{y_1, ..., y_N}\].
#' @param mu parameter \eqn{\mu}.
#' @param k parameter k.
#' @param theta parameter \eqn{\theta}.
#' @param rho parameter \eqn{\rho}.
#' @param h time unit.
#'
#' @return new value for \eqn{\sigma_v}, a scale value.
#' @export
#'
#' @examples
#' sigma_v_old = 0.01; prior_shape = 2; prior_rate = 1
#' v = c(0.25,0.2,0.1,0.15,0.25) # variance of diff(v) can't be 0!
#' y = rep(0.125,4)
#' mu = 0.125; k = 0.1; theta = 0.25; rho = 0
#' h = 1
#' rsigma_v(sigma_v_old,prior_shape,prior_rate,v,y,mu,k,theta,rho,h)
rsigma_v <- function(sigma_v_old,prior_shape,prior_rate,v,y,mu,k,theta,rho,h) {
  # check extreme case
  msg = 'variance of diff(v) = 0, thus proposal_rate = Inf'
  if(stats::var(diff(v)) == 0) stop(msg)
  #
  # propose a new sigma_v
  proposal_rate = 1/sqrt(stats::var(diff(v))/(mean(v)*h))
  sigma_v_new = stats::rgamma(1, shape = 2, rate= proposal_rate)
  ratio_proposal = exp(
    stats::dgamma(sigma_v_old, 2, rate=proposal_rate, log = TRUE) -
    stats::dgamma(sigma_v_new, 2, rate=proposal_rate, log = TRUE) )
  #
  sum_new = 0; sum_old = 0; N = length(y)
  for (n in 1:N) {
    # y_n <- v_{n-1}
    std_eps_y = y[n] - mu*h + v[n]*h/2
    # v_n <- v_{n-1}, y_n
    mu_new = sigma_v_new*rho*std_eps_y; var_new = sigma_v_new^2*(1-rho^2)*v[n]*h
    mu_old = sigma_v_old*rho*std_eps_y; var_old = sigma_v_old^2*(1-rho^2)*v[n]*h
    x_n = v[n+1] - v[n] - k * (theta - v[n]) * h
    #
    sum_new = sum_new + (x_n-mu_new)^2/(2*var_new)
    sum_old = sum_old + (x_n-mu_old)^2/(2*var_old)
  }
  ratio_lklhd = exp(-N * log(sigma_v_new/sigma_v_old) - sum_new + sum_old)
  #
  ratio_prior = exp(
    stats::dgamma(sigma_v_new, prior_shape, prior_rate, log=TRUE) -
    stats::dgamma(sigma_v_old, prior_shape, prior_rate, log=TRUE) )
  #
  ratio_post = ratio_lklhd * ratio_prior
  #
  ratio = ratio_post * ratio_proposal
  if(is.na(ratio)) {
    print(paste0('sigma_v_old: ',sigma_v_old,'; sigma_v_new:',sigma_v_new))
    print(paste0('ratio_lklhd: ',ratio_lklhd))
    print(paste0('ratio_prior: ',ratio_prior))
    print(paste0('ratio_proposal: ',ratio_proposal))
    stop('NaN produced in rsigma_v()!')
  }
  #
  p_accept = min(ratio,1)
  U = stats::runif(1)
  if (U <= p_accept) {
    sigma_v = sigma_v_new
  } else {
    sigma_v = sigma_v_old
  }
  #
  return(sigma_v)
}
