#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rmu
// [[Rcpp::export]]
double crmu(double prior_mu, double prior_var, NumericVector v, NumericVector y,
            double k, double theta, double sigma_v, double rho, double h) {
  int N = y.size();
  double std_eps_v, x_n, lklhd_var, post_mu, post_var;
  for(int n=0; n < N; n++) {
    // v_n <- v_{n-1}
    std_eps_v = (v[n+1] - v[n] - k * (theta - v[n]) * h )/sigma_v;
    // y_n <- v_{n-1}, v_n
    x_n = (y[n] + v[n]*h/2 - rho*std_eps_v ) / h;
    lklhd_var = (1- rho*rho) * v[n] / h;
    //
    // update parameters, posterior mu and variance
    post_var = 1/(1/prior_var + 1/lklhd_var);
    post_mu = post_var * (prior_mu/prior_var + x_n/lklhd_var);
    //
    prior_mu = post_mu;
    prior_var = post_var;
  }
  if (post_var < 0) {
    Rcerr << "crmu: post_var < 0!" << "\n";
  }
  double mu = R::rnorm(post_mu, std::sqrt(post_var));
  return mu;
}
