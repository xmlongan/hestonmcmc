#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
using namespace Rcpp;

double cor(NumericVector v1, NumericVector v2) {
  int N = v1.size();
  double cov =0;
  double mean1 = Rcpp::mean(v1);
  double mean2 = Rcpp::mean(v2);
  for(int n=0; n < N; n++) {
    cov += (v1[n]-mean1) * (v2[n]-mean2);
  }
  cov = cov / (N-1);
  double rho0 = cov / (Rcpp::sd(v1) * Rcpp::sd(v2));
  return rho0;
}

double cget_rho0(NumericVector v, NumericVector y, double mu, double k,
                double theta, double sigma_v, double h) {
  int N = y.size();
  NumericVector eps_v (N);
  NumericVector eps_y (N);
  double rho0;
  for(int n=0; n < N; n++) {
    eps_y[n] = (y[n] - mu * h + v[n] * h/2) / std::sqrt(v[n] * h);
    eps_v[n] = (v[n+1]-v[n] - k*(theta-v[n])*h)/(sigma_v * std::sqrt(v[n] *h));
  }
  if (Rcpp::sd(eps_y) == 0 || Rcpp::sd(eps_v) == 0) {
    rho0 = 0;
    Rcout << "crrho -> crrho_proposal -> cget_rho0: ";
    Rcout << "sd(eps_y) or sd(eps_v) equals 0!\n";
  } else {
    rho0 = cor(eps_y, eps_v);
  }
  //
  if (NumericVector::is_na(rho0)) {
    Rcpp::Rcerr << "cget_rho0: rho0 is NaN!\n";
  }
  return rho0;
}

double cdrho_proposal(double x, double rho0) {
  double d;
  if (rho0 != -1 && rho0 != 1) {
    if (x < rho0) {
      d = (x + 1) / (rho0 + 1);
    } else {
      d = (1 - x) / (1 - rho0);
    }
  } else {
    d = 1/2 + rho0 * x / 2;
  }
  return d;
}

double crrho_proposal(double rho0) {
  double rho;
  double u = R::runif(0,1);
  if (u < (rho0 + 1)/2) {
    rho = std::sqrt(2 * u * (rho0+1)) - 1;
  } else {
    rho = 1 - std::sqrt(2 * (1-u) * (1-rho0));
  }
  return rho;
}


//' @rdname rrho
// [[Rcpp::export]]
double crrho(double rho_old, NumericVector v, NumericVector y, double mu,
             double k, double theta, double sigma_v, double h) {
  double rho0, rho_new, ratio_proposal;
  double sum_old = 0, sum_new = 0, std_eps_v, x_n, ratio_post, ratio, p_accept;
  int N = y.size();
  // compute the sample correlation between the Brownian increments
  rho0 = cget_rho0(v, y, mu, k, theta, sigma_v, h);
  // propose a new rho
  rho_new = crrho_proposal(rho0);
  ratio_proposal = cdrho_proposal(rho_old, rho0) / 
                   cdrho_proposal(rho_new, rho0);
  //
  for(int n=0; n < N; n++) {
    // y_n <- v_{n-1},v_n
    x_n = y[n] - mu*h + v[n]*h/2;
    // v_n <- v_{n-1}
    std_eps_v = (v[n+1] - v[n] - k * (theta - v[n]) * h ) / sigma_v;
    sum_old += std::pow(x_n - rho_old * std_eps_v, 2) /
               (2 * (1 - rho_old * rho_old) * v[n] * h);
    sum_new += std::pow(x_n - rho_new * std_eps_v, 2) /
               (2 * (1 - rho_new * rho_new) * v[n] * h);
  }
  ratio_post = std::exp(-std::log((1-rho_new*rho_new)/(1-rho_old*rho_old))*N/2
               - sum_new + sum_old);
  //
  ratio = ratio_post * ratio_proposal;
  //
  if (NumericVector::is_na(ratio) ) {
    Rcerr << "crrho: ratio is NaN!\n";
  }
  //
  p_accept = (ratio < 1) ? ratio : 1;
  double rho = rho_old;
  if (R::runif(0,1) <= p_accept) {
    rho = rho_new;
  }
  return rho;
}
