#include <Rcpp.h>
#include <cmath>

#include "crmu.h"
#include "crk.h"
#include "crtheta.h"
#include "crsigma_v.h"
#include "crrho.h"

#include "crv0.h"
#include "crvi.h"
#include "crvN.h"

#include "cinitialize_values.h"

using namespace Rcpp;

//' @rdname mcmc
// [[Rcpp::export]]
NumericVector cmcmc(NumericVector y, int g = 5000, int G = 10000, double h = 1,
                    double echo = 0) {
  // g: warm-up samples
  // G: total samples
  NumericVector parameters (5);
  NumericVector v = cinitialize_values(y,parameters,1);
  //
  double mu, k, theta, sigma_v, rho;
  mu = parameters[0];
  k = parameters[1];
  theta = parameters[2];
  sigma_v = parameters[3];
  rho = parameters[4];
  //
  double prior_shape_v0, prior_rate_v0;
  prior_shape_v0 = 2*k*theta/std::pow(sigma_v, 2);
  prior_rate_v0  = 2*k/std::pow(sigma_v, 2);
  //
  int N = y.size();
  //
  double cum_mu=0, cum_k=0, cum_theta=0, cum_sigma_v=0, cum_rho=0;
  //
  double prior_mu, prior_var, prior_shape, prior_rate;
  for (int i=0; i < G; i++) {
    // update parameters
    //
    prior_mu  = 0;
    prior_var = 0.1*0.1;
    mu = crmu(prior_mu, prior_var, v, y, k, theta, sigma_v, rho, h);
    //
    prior_mu  = 0.01;
    prior_var = 1*1;
    k = crk(prior_mu, prior_var, v, y, mu, theta, sigma_v, rho, h);
    //
    prior_mu  = 0.01;
    prior_var = 1*1;
    theta = crtheta(prior_mu, prior_var, v, y, mu, k, sigma_v, rho, h);
    //
    prior_shape = 2;
    prior_rate  = 1;
    sigma_v = crsigma_v(sigma_v, prior_shape, prior_rate,v,y,mu,k,theta,rho,h);
    //
    rho = crrho(rho, v, y, mu, k, theta, sigma_v, h);
    //
    // update v0
    v[0] = crv0(v[0], prior_shape_v0, prior_rate_v0, v[1], y[0], mu, k, theta,
                sigma_v, rho, h);
    // update v1:v_{N-1}
    for (int n=1; n < N; n++) {
      v[n] = crvi(v[n], v[n-1], v[n+1], y[n-1], y[n], mu,k,theta,sigma_v,rho,h);
    }
    // update vN
    v[N] = crvN(v[N], v[N-1], y[N-1], mu, k, theta, sigma_v, rho, h);
    //
    // cumulate parameters
    if (i > g-1) {
      cum_mu += mu;
      cum_k  += k;
      cum_theta += theta;
      cum_sigma_v += sigma_v;
      cum_rho += rho;
    }
    //
    if (echo == 1 && (i < 10 || i > G-11)) {
      Rcout << i+1 << "th:\nparameters: " << mu << "," << k << ",";
      Rcout << theta << "," << sigma_v << "," << rho << "\n";
      Rcout << "first 10 of v: " << v[0] << "," << v[1] << "," << v[2] << ",";
      Rcout << v[3] << "," << v[4] << "," << v[5] << "," << v[6] << ",";
      Rcout << v[7] << "," << v[8] << "," << v[9] << "\n";
    }
  }
  // estimate
  int n = G - g;
  NumericVector estimate={cum_mu/n,cum_k/n,cum_theta/n,cum_sigma_v/n,cum_rho/n};
  return estimate;
}
