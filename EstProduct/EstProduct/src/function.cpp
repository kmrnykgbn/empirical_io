#include <stdio.h>
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]
/** moment_0P_2nd <- function(alpha, beta_0, beta_k, df, df_1st) {
  moment = (1/(max(df$t)*max(df$j,na.rm =TRUE)))*sum(df_1st$y_error_tilde-beta_0-beta_k*df$k
                                                       -alpha*(df_1st$phi_t_1-beta_0-beta_k*df$k))
  *lag()
  return(moment)
} **/
// [[Rcpp::export]]
NumericVector moment_OP_2nd_rcpp(int T, int J, double alpha, double beta_0, double beta_k,
                                 Eigen::VectorXd y_error_tilde, Eigen::VectorXd phi_t_1, 
                                 Eigen::MatrixXd matA) {
  NumericVector r = {0.0, 0.0, 0.0};
  
  
  for (int j = 0; j < J; j++) {
    for (int t = 0; t < T; t++) {
      const int i = j*T + t;
      const int K = matA.cols();
      
      double k = matA(i , 0);
      double lag_k = matA(i , 1);
      
      for (int k = 0; k < K; k++) {
        double resVal = (y_error_tilde(i) -beta_0 - beta_k * k 
          -alpha * (phi_t_1(i) -beta_0 - beta_k * lag_k))
          * matA(i , k);
        r[k] += resVal;
      }
    }
  }
  r = r/(J*T);
  return r;
}

  