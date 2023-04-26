#include <Rcpp.h>
using namespace Rcpp;

double crtheta(double prior_mu, double prior_var, NumericVector v,
               NumericVector y, double mu, double k, double sigma_v, double rho,
               double h);
