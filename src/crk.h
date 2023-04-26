#include <Rcpp.h>
using namespace Rcpp;

double crk(double prior_mu, double prior_var, NumericVector v, NumericVector y,
           double mu, double theta, double sigma_v, double rho, double h);
