#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rk
// [[Rcpp::export]]
double crk(double prior_mu, double prior_var, NumericVector v, NumericVector y,
           double mu, double theta, double sigma_v, double rho, double h) {
  int N = y.size();
  double std_eps_y, x_n, lklhd_var, post_mu, post_var;
  for(int n=0; n < N; n++) {
    // y_n <- v_{n-1}
    std_eps_y = y[n] - mu*h + v[n]*h/2;
    //
    // v_n <- v_{n-1},y_n
    x_n = v[n+1] - v[n] - sigma_v * rho * std_eps_y;
    lklhd_var = sigma_v * sigma_v * (1-rho*rho) * v[n] * h;
    //
    // update parameters, posterior mu and variance
    post_var = 1/(1/prior_var + std::pow((theta-v[n])*h, 2)/lklhd_var);
    post_mu = post_var * (prior_mu/prior_var + x_n*(theta-v[n])*h/lklhd_var);
    //
    prior_mu = post_mu;
    prior_var = post_var;
  }
  if (post_var < 0) {
    Rcerr << "crk: post_var < 0!\n";
  }
  double k = R::rnorm(post_mu,std::sqrt(post_var));
  return k;
}
