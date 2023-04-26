#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname initialize_values
// [[Rcpp::export]]
NumericVector cinitialize_values(NumericVector y, NumericVector parameters,
                                 double h) {
  double std_eps_y, mu_v, sd_v;
  // initialize parameter values
  // mu=0, k=0.01, theta=0.1, sigma_v=0.01, rho=0
  parameters[0] = 0;
  parameters[1] = 0.01;
  parameters[2] = 0.1;
  parameters[3] = 0.01;
  parameters[4] = 0;
  double mu=0, k=0.01, theta=0.1, sigma_v=0.01, rho=0;
  //
  int N = y.size();
  NumericVector v (N+1);
  // let v_0 be the long-run mean of volatility
  v[0] = theta;
  //
  for (int n=1; n < N+1; n++) {
    // y_n <- v_{n-1}
    std_eps_y = y[n-1] - mu*h + v[n-1]*h/2;
    // v_n <- v_{n-1},y_n
    mu_v = v[n-1] + k * (theta - v[n-1]) * h + sigma_v * rho * std_eps_y;
    sd_v = sigma_v * std::sqrt((1-rho*rho) * v[n-1] * h);
    v[n] = R::rnorm(mu_v, sd_v);
    if (v[n] <= 0) {
      v[n] = 0.00001;
    }
  }
  return v;
}
