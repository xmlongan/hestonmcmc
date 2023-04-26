#' Estimate the parameters through MCMC
#'
#' @description 
#' Given \eqn{[y_1,\cdots, y_N]}, to estimate parameters of Heston Stochastic
#' Volatility Model using MCMC.
#' * `cmcmc()` implemented in C++, returns only a vector of the estimated
#' parameters.
#' * `mcmc()` implemented in R.
#'
#' @details 
#' Priors for and updating functions of these parameters:
#' |prior|updating function|
#' |:---|:---|
#' |\eqn{\mu \sim \mathcal{N}(0,0.1^2) }|  [rmu()]|
#' |\eqn{ k \sim \mathcal{N}(0.01,1^2) }|  [rk()]|
#' |\eqn{\theta \sim \mathcal{N}(0.01,1^2)}|  [rtheta()]|
#' |\eqn{\sigma_v \sim Gamma(shape=2,rate=1)}|  [rsigma_v()]|
#' |\eqn{\rho \sim Uniform(-1,1)}|  [rrho()]|
#'
#' Prior for the starting volatility: an uniform distribution over a large 
#' range. (Depreciated 
#' \eqn{v_0 \sim Gamma(shape=2k\theta/\sigma_v^2, rate=2k/\sigma_v^2)}.)
#'
#' Updating functions of volatility:
#' * \eqn{v_0}: [rv0()]
#' * \eqn{v_n}: [rvi()] for \eqn{1\le n < N}
#' * \eqn{v_N}: [rvN()]
#'
#' @param y vector of returns \[\eqn{y_1, ..., y_N}\].
#' @param g number of warm-up iteration, defaults to 5,000.
#' @param G number of total iteration, defaults to 10,000.
#' @param h time unit, defaults to 1.
#' @param echo TRUE or FALSE, whether or not echo the parameters and first 10
#' volatilities in the first and last ten iterations, defaults to FALSE.
#'
#' @return a list of two elements:
#' * a vector of the estimated parameters, \eqn{(\mu,k,\theta,\sigma_v,\rho)},
#' * a matrix record of volatility in which each row represents an
#' iteration result, i.e., with dimension of \eqn{G\times(N+1)}.
#' @export
#'
#' @examples
#' y = rep(0.125,20)
#' mcmc(y)
mcmc <- function(y,g=5000,G=10000,h=1,echo=FALSE) {
  # g: warm-up samples
  # G: total samples
  init_values = initialize_values(y, parameters=c(0,0.01,0.1,0.01,0), h=1)
  #
  mu      = init_values$parameters[1]
  k       = init_values$parameters[2]
  theta   = init_values$parameters[3]
  sigma_v = init_values$parameters[4]
  rho     = init_values$parameters[5]
  #
  prior_shape_v0 = 2*k*theta/sigma_v^2
  prior_rate_v0 = 2*k/sigma_v^2
  #
  v = init_values$v
  #
  N = length(y)
  #
  record_mu = rep(0,G)
  record_k  = rep(0,G)
  record_theta = rep(0,G)
  record_sigma_v = rep(0,G)
  record_rho = rep(0,G)
  record_v = matrix(rep(0,G*(N+1)), nrow=G, ncol=N+1)
  #
  for (i in 1:G) {
    # update parameters
    #
    prior_mu = 0; prior_var = 0.1^2
    mu = rmu(prior_mu, prior_var, v, y, k, theta, sigma_v, rho, h)
    #
    prior_mu = 0.01; prior_var = 1^2
    k = rk(prior_mu, prior_var, v, y, mu, theta, sigma_v, rho, h)
    #
    prior_mu = 0.01; prior_var = 1^2
    theta = rtheta(prior_mu, prior_var, v, y, mu, k, sigma_v, rho, h)
    #
    prior_shape = 2; prior_rate = 1
    sigma_v = rsigma_v(sigma_v, prior_shape, prior_rate, v,y,mu,k,theta,rho,h)
    #
    rho = rrho(rho, v, y, mu, k, theta, sigma_v, h)
    #
    # update v0
    v[1] = rv0(v[1], prior_shape_v0, prior_rate_v0, v[2], y[1],
               mu, k, theta, sigma_v, rho, h)
    record_v[i,1] = v[1]
    # update  v1:v_{N-1}
    for (n in 1:(N-1)) {
      v[n+1] = rvi(v[n+1], v[n], v[n+2], y[n], y[n+1],
                   mu, k, theta, sigma_v, rho, h)
      record_v[i,n+1] = v[n+1]
    }
    # update vN
    v[N+1] = rvN(v[N+1], v[N], y[N], mu, k, theta, sigma_v, rho, h)
    record_v[i,N+1] = v[N+1]
    #
    # record
    record_mu[i] = mu
    record_k[i] = k
    record_theta[i] = theta
    record_sigma_v[i] = sigma_v
    record_rho[i] = rho
    #
    if (echo==TRUE && (i<11 || i>G-10)) {
      cat(sprintf('\n%dth: \nparameters: ',i))
      cat(signif(c(mu,k,theta,sigma_v,rho), digits=4), sep=',')
      cat('\nfrist ten of v: ')
      cat(signif(v[1:10], digits = 3), sep = ',')
    }
  }
  #
  # estimate
  mu = mean(record_mu[g:G])
  k  = mean(record_k[g:G])
  theta = mean(record_theta[g:G])
  sigma_v = mean(record_sigma_v[g:G])
  rho = mean(record_rho[g:G])
  #
  return(list(parameters = c(mu,k,theta,sigma_v,rho),
              record_v = record_v[seq(G-49,G), seq(N-9,N)]))
}
