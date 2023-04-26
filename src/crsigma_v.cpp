#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rsigma_v
// [[Rcpp::export]]
double crsigma_v(double sigma_v_old, double prior_shape, double prior_rate,
                 NumericVector v, NumericVector y, double mu, double k,
                 double theta, double rho, double h) {
  int N = y.size();
  double sigma_v_new, ratio_proposal, sum_new, sum_old;
  double mu_new, mu_old, var_new, var_old;
  double std_eps_y, x_n, ratio_lklhd, ratio_prior, ratio_post, ratio, p_accept;
  //
  double S = Rcpp::var(Rcpp::diff(v));
  if(S == 0) {
    Rcerr << "crsigma_v: var(diff(v)) = 0, therefore scale of proposal is 0!";
    Rcerr << "\n";
  }
  // propose a new sigma_v
  double proposal_scale = std::sqrt(S / (Rcpp::mean(v)*h) );
  sigma_v_new = R::rgamma(2, proposal_scale);// shape, scale
  ratio_proposal = std::exp(R::dgamma(sigma_v_old,2,proposal_scale,1) -
                            R::dgamma(sigma_v_new,2,proposal_scale,1));
  //
  sum_new = 0;
  sum_old = 0;
  for(int n=0; n < N; n++) {
    // y_n <- v_{n-1}
    std_eps_y = y[n] - mu*h + v[n]*h/2;
    // v_n <- v_{n-1},y_n
    x_n = v[n+1] - v[n] - k * (theta - v[n]) * h;
    mu_new = sigma_v_new * rho * std_eps_y;
    mu_old = sigma_v_old * rho * std_eps_y;
    var_new = std::pow(sigma_v_new,2) * (1-rho*rho) * v[n] * h;
    var_old = std::pow(sigma_v_old,2) * (1-rho*rho) * v[n] * h;
    //
    sum_new += std::pow(x_n-mu_new, 2) / (2*var_new);
    sum_old += std::pow(x_n-mu_old, 2) / (2*var_old);
  }
  ratio_lklhd = std::exp( -N * std::log(sigma_v_new/sigma_v_old)
                         - sum_new + sum_old);
  //
  ratio_prior = std::exp(R::dgamma(sigma_v_new, prior_shape, 1/prior_rate, 1) -
                         R::dgamma(sigma_v_old, prior_shape, 1/prior_rate, 1));
  //
  ratio_post = ratio_lklhd * ratio_prior;
  //
  ratio = ratio_post * ratio_proposal;
  //
  if(NumericVector::is_na(ratio) ) {
    Rcerr << "crsigma_v: ratio is NaN!\n";
  }
  //
  p_accept = (ratio < 1) ? ratio : 1;
  double sigma_v = sigma_v_old;
  if(R::runif(0,1) <= p_accept) {
    sigma_v = sigma_v_new;
  }
  return sigma_v;
}
