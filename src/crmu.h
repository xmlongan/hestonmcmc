#include <Rcpp.h>
using namespace Rcpp;

double crmu(double prior_mu, double prior_var, NumericVector v, NumericVector y,
            double k, double theta, double sigma_v, double rho, double h);
