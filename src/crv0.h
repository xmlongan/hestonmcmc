#include <Rcpp.h>
using namespace Rcpp;

double crv0(double v0_old, double prior_shape, double prior_rate, double v1,
            double y1, double mu, double k, double theta, double sigma_v,
            double rho, double h);
