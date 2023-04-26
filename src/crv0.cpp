#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rv0
// [[Rcpp::export]]
double crv0(double v0_old, double prior_shape, double prior_rate, double v1,
            double y1, double mu, double k, double theta, double sigma_v,
            double rho, double h) {
  double v0_new, ratio_proposal, ratio_prior, ratio1, m_new, m_old, v_new,v_old;
  double ratio2, ratio_post, ratio, p_accept;
  // propose new v0
  v0_new = R::rgamma(2,v0_old);// shape,scale
  ratio_proposal = std::exp( R::dgamma(v0_old, 2, v0_new, 1) -
                             R::dgamma(v0_new, 2, v0_old, 1));
  ratio_prior = 1;// assuming uniform distribution for the prior
  //
  // y_1 <- v_0
  ratio1 = std::exp( -std::pow(y1 - mu*h + v0_new*h/2, 2) / (2*v0_new*h) +
                      std::pow(y1 - mu*h + v0_old*h/2, 2) / (2*v0_old*h) );
  //
  // v_1 <- v_0,y_1
  m_new = v0_new + k*(theta-v0_new)*h + sigma_v*rho*(y1 - mu*h + v0_new*h/2);
  m_old = v0_old + k*(theta-v0_old)*h + sigma_v*rho*(y1 - mu*h + v0_old*h/2);
  v_new = sigma_v * sigma_v * (1-rho*rho) * v0_new * h;
  v_old = sigma_v * sigma_v * (1-rho*rho) * v0_old * h;
  ratio2 = std::exp( -std::pow(v1-m_new, 2) / (2*v_new) +
                      std::pow(v1-m_old, 2) / (2*v_old));
  //
  ratio_post = (v0_old/v0_new) * ratio2 * ratio1 * ratio_prior;
  //
  ratio = ratio_post * ratio_proposal;
  if (NumericVector::is_na(ratio)) {
    Rcerr << "crv0: ratio is NaN!\n";
  }
  //
  p_accept = (ratio < 1) ? ratio : 1;
  double v0 = v0_old;
  if (R::runif(0,1) <= p_accept) {
    v0 = v0_new;
  }
  return v0;
}
