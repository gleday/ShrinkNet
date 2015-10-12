#include <RcppArmadillo.h>
#include <R.h>

// [[Rcpp::export]]
arma::colvec mydigamma(arma::colvec vec){
  arma::colvec out(vec.n_elem);
  for(int k = 0; k<vec.n_elem; k++){
    out(k) = R::digamma(vec(k));
  }
  return out;
}
