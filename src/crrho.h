#include <Rcpp.h>
using namespace Rcpp;

double cor(NumericVector v1, NumericVector v2);
double get_rho0(NumericVector v, NumericVector y, double mu, double k,
                double theta, double sigma_v, double h);
double drho_proposal(double x, double rho0);
double rrho_proposal(double rho0);
double crrho(double rho_old, NumericVector v, NumericVector y, double mu,
             double k, double theta, double sigma_v, double h);
