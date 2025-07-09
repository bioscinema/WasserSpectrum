#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List compute_beta1_se(NumericMatrix Phi,
                      NumericVector coef,
                      NumericMatrix cov_matrix,
                      IntegerVector sel_idx) {
  int m = Phi.nrow();
  int K = Phi.ncol();
  NumericVector beta1(m);
  NumericVector se1(m);
  
  for (int j = 0; j < m; ++j) {
    double b = 0.0;
    for (int k = 0; k < K; ++k) {
      b += coef[sel_idx[k]] * Phi(j, k);
    }
    beta1[j] = b;
    
    // compute variance
    double var = 0.0;
    for (int k1 = 0; k1 < K; ++k1) {
      for (int k2 = 0; k2 < K; ++k2) {
        var += Phi(j, k1) * cov_matrix(sel_idx[k1], sel_idx[k2]) * Phi(j, k2);
      }
    }
    se1[j] = sqrt(var);
  }
  
  return List::create(
    Named("beta1") = beta1,
    Named("se1") = se1
  );
}

