#' Initialize the Latent Volatility and Parameters
#'
#' @description 
#' Initialize the latent volatility and parameters with a good guess. Making
#' sure the initialized volatility are positive.
#' * `cinitialize_values()` implemented in C++, through package Rcpp, returns a
#' vector of volatility, the outside `parameters` changes along with the 
#' initialization inside the function.
#' * `initialize_values()` implemented in R.
#'
#' @details 
#' Parameters are first initialized. The volatility starting
#' point \eqn{v_0} is set to the long-run mean of volatility, \eqn{\theta}.
#' Then we use \eqn{P(v_n|v_{n-1},y_n)} to sample \eqn{v_n} iteratively.
#' However, we may end up some non-positive volatility because of using the
#' discretized Heston SV model, i.e.,
#' \deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}\epsilon_n^y,}
#' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}\epsilon_n^v.}
#' The conditional distribution
#' \deqn{P(v_n|v_{n-1},y_n)\sim\mathcal{N}
#'       (v_{n-1}+k(\theta-v_{n-1})h + \sigma_v\rho\sqrt{v_{n-1}h}\epsilon_n^y,
#'        \sigma_v^2(1-\rho^2)v_{n-1}h)}
#' where \deqn{\sqrt{v_{n-1}h}\epsilon_n^y = y_n-\mu h + \frac{1}{2}v_{n-1}h.}
#' If that does happen, we make a change as
#' \deqn{v_n =
#'       \begin{cases}
#'        v_n|v_{n-1},y_n & \text{if } v_n|v_{n-1},y_n > 0,\\
#'        0.00001         & \text{if } v_n|v_{n-1},y_n \le 0.
#'       \end{cases}}
#'
#' @param y vector of returns \[\eqn{y_1, ..., y_N}\].
#' @param parameters vector of parameters \eqn{(\mu,k,\theta,\sigma_v,\rho)}.
#' @param h time unit.
#'
#' @return a list with two elements, one vector of all parameters, and another
#' vector of the latent volatility.
#' @export
#'
#' @examples
#' y = rep(0.125,4)
#' parameters_v = initialize_values(y)
initialize_values <- function(y,parameters=c(0,0.01,0.1,0.01,0),h=1) {
  # initial parameter values
  # mu=0; k=0.01; theta=0.1; sigma_v=0.01; rho=0
  mu = parameters[1]; k = parameters[2]; theta = parameters[3]
  sigma_v = parameters[4]; rho = parameters[5]

  N = length(y)
  # initial v0:vN
  v = rep(0,N)
  # v0 ~ Gamma(shape,rate)
  # v[1] = rgamma(1, shape=2*k*theta/sigma_v^2, rate=2*k/sigma_v^2)
  # however, let it be the long-run mean of volatility
  v[1] = theta
  #
  for (n in 1:N) {
    # y_n <- v_{n-1}
    std_eps_y = y[n] - mu*h + v[n]*h/2
    # v_n <- v_{n-1},y_n
    mu_v = v[n] + k*(theta-v[n])*h + sigma_v * rho * std_eps_y
    sd_v = sigma_v * sqrt((1-rho^2) * v[n]*h)
    v[n+1] = stats::rnorm(1, mean = mu_v, sd = sd_v)
    if (v[n+1] <= 0) {
      v[n+1] = 0.00001
    }
  }
  return(list(parameters=c(mu,k,theta,sigma_v,rho),v=v))
}
