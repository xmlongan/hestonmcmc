#include <Rcpp.h>
#include <Rmath.h>
#include <cmath>
#include <string>
#include <fstream>
using namespace Rcpp;

//' Generate Heston SV Sample Trajectory
//' 
//' @description
//' Generate Heston Stochastic Volatility sample trajectory using Euler 
//' Approximation.
//' 
//' @details
//' Heston Stochastic Volatility:
//' \deqn{y_n=\mu h -\frac{1}{2}v_{n-1}h + \sqrt{v_{n-1}h}(\rho\epsilon_n^v
//'       +\sqrt{1-\rho^2}\epsilon_n),}
//' \deqn{v_n-v_{n-1}=k(\theta-v_{n-1})h + \sigma_v\sqrt{v_{n-1}h}\epsilon_n^v.}
//' 
//' @param v_0 the initial volatility \eqn{v_0}.
//' @param n_segment number of further segments for each interval.
//' @param par a vector of true values for the parameters
//' \eqn{\mu,k,\theta,\sigma_v,\rho}.
//' @param N the number of sample returns, \eqn{y_1,\cdots,y_N}.
//' @param h time unit.
//' 
//' @return a vector of returns, \eqn{(y_1,\cdots,y_N)}.
//' @export
//' 
//' @examples
//' S0 = c(0.125,0.1,0.25,0.1,-0.7)
//' y_series_0 = crHeston(v_0=S0[3], n_segment=10, par=S0, N=1000, h=1)
//' 
//' S1 = c(0.4,0.1,0.25,0.1,-0.7)
//' y_series_1 = crHeston(v_0=S1[3], n_segment=10, par=S1, N=1000, h=1)
//' 
//' S2 = c(0.125,0.03,0.25,0.1,-0.7)
//' y_series_2 = crHeston(v_0=S2[3], n_segment=10, par=S2, N=1000, h=1)
//' 
//' S3 = c(0.125,0.1,0.5,0.1,-0.7)
//' y_series_3 = crHeston(v_0=S3[3], n_segment=10, par=S3, N=1000, h=1)
//' 
//' S4 = c(0.125,0.1,0.25,0.2,-0.7)
//' y_series_4 = crHeston(v_0=S4[3], n_segment=10, par=S4, N=1000, h=1)
//' 
//' S5 = c(0.125,0.1,0.25,0.1,-0.3)
//' y_series_5 = crHeston(v_0=S5[3], n_segment=10, par=S5, N=1000, h=1)
// [[Rcpp::export]]
NumericVector crHeston(double v_0,double n_segment,NumericVector par,int N,
                       double h=1){
  double mu,k,theta,sigma,rho;
  mu = par[0]; k = par[1]; theta = par[2]; sigma = par[3]; rho = par[4];
  NumericVector v_series (N+1); NumericVector y_series (N);
  v_series[0] = v_0;
  double I, I2, IV, eps, eps2, dlt, v, y;
  dlt = h/n_segment; // dlt: delta
  
  for (int i = 0; i < N; i++) {
    I = 0; I2 = 0; IV = 0;
    for (int j = 0; j < n_segment; j++) {
      eps  = std::sqrt(v_0) * R::rnorm(0,std::sqrt(dlt));
      eps2 = std::sqrt(v_0) * R::rnorm(0,std::sqrt(dlt));
      v    = v_0 + dlt*k*(theta-v_0) + sigma*eps;
      if (v < 0) {
        v = 0;
        eps = (-v_0-dlt*k*(theta-v_0))/sigma;
      }
      I += eps; I2 += eps2; IV += v_0*dlt;
      v_0 = v;
    }
    v_series[i+1] = v;
    y = mu*h - IV/2 + rho*I + std::sqrt(1-rho*rho)*I2;
    y_series[i] = y;
  }
  return y_series;
}


void write_csv(std::string tofile,NumericVector y_series) {
  std::ofstream dfile(tofile);
  for (int i = 0; i < y_series.size(); i++) {
    dfile << y_series[i] << "\n";
  }
  dfile.close();
}


//' Generate and Write to Files Large Number of Sample Trajectories
//' 
//' @description
//' Generate and Write to Files Large Number of Sample Trajectories. The files
//' are named according to 'par_name-N-ith.csv' where N denotes sample length.
//' 
//' @param par_name parameter setting name, such as 'par0'.
//' @param par a vector of true values for the parameters
//' \eqn{\mu,k,\theta,\sigma_v,\rho}.
//' @param N the number of sample returns, \eqn{y_1,\cdots,y_N}.
//' @param N_rep the number of Sample Trajectories.
//' @param h time unit, defaults to 1.
//' @param n_segment number of further segments for each interval, 
//' defaults to 10.
//' 
//' @return no return.
//' 
//' @examples
//' # par0 = c(0.125,0.1,0.25,0.1,-0.7)
//' # gen_data('par0', par0, N=100000, N_rep=0, h=1, n_segment=10)
// [[Rcpp::export]]
void gen_data(std::string par_name,NumericVector par,int N,int N_rep,
              double h=1,double n_segment=10) {
  std::string tofile;
  std::string prefix;
  prefix = par_name + "-" + std::to_string(N);
  double theta = par[2];
  NumericVector y_series;
  for (int i = 0; i < N_rep; i++) {
    tofile = prefix + "-" + std::to_string(i+1) + "th.csv";
    y_series = crHeston(theta,n_segment,par,N,h);
    write_csv(tofile,y_series);
  }
}
