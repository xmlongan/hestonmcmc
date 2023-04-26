#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

//' @rdname rvN
// [[Rcpp::export]]
double crvN(double vN_old, double vNm1, double yN, double mu, double k,
            double theta, double sigma_v, double rho, double h) {
  double vN_new, ratio_proposal, muN, varN, ratio_post, ratio, p_accept, vN;
  // propose new vN
  vN_new = R::rgamma(2, vN_old);// shape, scale
  ratio_proposal = std::exp( R::dgamma(vN_old, 2, vN_new, 1) -
                             R::dgamma(vN_new, 2, vN_old, 1));
  //
  // v_N <- v_{N-1},y_N
  muN  = vNm1 + k*(theta-vNm1)*h + sigma_v*rho*(yN - mu*h + vNm1*h/2);
  varN = sigma_v * sigma_v * (1-rho*rho) * vNm1 * h;
  //
  ratio_post = std::exp( -std::pow(vN_new-muN, 2) / (2*varN) +
                          std::pow(vN_old-muN, 2) / (2*varN) );
  //
  ratio = ratio_post * ratio_proposal;
  if (NumericVector::is_na(ratio)) {
    Rcerr << "crvN: ratio is NaN!\n";
  }
  //
  p_accept = (ratio < 1) ? ratio : 1;
  vN = vN_old;
  if (R::runif(0,1) <= p_accept) {
    vN = vN_new;
  }
  return vN;
}
