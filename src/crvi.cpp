#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rvi
// [[Rcpp::export]]
double crvi(double vi_old, double vim1, double vip1, double yi, double yip1,
            double mu, double k, double theta, double sigma_v, double rho,
            double h) {
  double vi_new, ratio_proposal, mu1, var1, ratio1, ratio2, mu3_new, mu3_old;
  double var_new, var_old, ratio3, ratio_post, ratio, p_accept, vi;
  // propose new vi
  vi_new = R::rgamma(2,vi_old);// shape, scale
  //
  ratio_proposal = std::exp( R::dgamma(vi_old, 2, vi_new, 1) -
                             R::dgamma(vi_new, 2, vi_old, 1));
  //
  // v_i <- v_{i-1},y_i
  mu1  = vim1 + k*(theta-vim1)*h + sigma_v*rho*(yi - mu*h + vim1*h/2);
  var1 = sigma_v * sigma_v * (1-rho*rho) * vim1 * h;
  ratio1 = std::exp( - std::pow(vi_new - mu1, 2) / (2*var1) +
                       std::pow(vi_old - mu1, 2) / (2*var1) );
  //
  // y_{i+1} <- v_i
  ratio2 = std::exp( - std::pow(yip1 - mu*h + vi_new*h/2, 2) / (2*vi_new*h) +
                       std::pow(yip1 - mu*h + vi_old*h/2, 2) / (2*vi_old*h) );
  //
  // v_{i+1} <- v_i,y_{i+1}
  mu3_new = vi_new + k*(theta-vi_new)*h + sigma_v*rho*(yip1 - mu*h+vi_new*h/2);
  mu3_old = vi_old + k*(theta-vi_old)*h + sigma_v*rho*(yip1 - mu*h+vi_old*h/2);
  var_new = sigma_v * sigma_v * (1-rho*rho) * vi_new * h;
  var_old = sigma_v * sigma_v * (1-rho*rho) * vi_old * h;
  ratio3 = std::exp( - std::pow(vip1-mu3_new, 2) / (2*var_new) +
                       std::pow(vip1-mu3_old, 2) / (2*var_old) );
  //
  ratio_post = (vi_old/vi_new) * ratio3 * ratio2 * ratio1;
  //
  ratio = ratio_post * ratio_proposal;
  if (NumericVector::is_na(ratio)) {
    Rcerr << "crvi: ratio is NaN!\n";
  }
  //
  p_accept = (ratio < 1) ? ratio : 1;
  vi = vi_old;
  if (R::runif(0,1) <= p_accept) {
    vi = vi_new;
  }
  return vi;
}
