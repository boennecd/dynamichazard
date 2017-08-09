#include "sample_funcs.h"
#include <RcppArmadilloExtensions/sample.h>

arma::uvec sample_indices(arma::vec probs){
  return sample_indices(probs.n_elem, probs);
}

arma::uvec sample_indices(int size, arma::vec probs){
  arma::uvec idx = arma::linspace<arma::uvec>(0, probs.n_elem - 1, probs.n_elem);

  return Rcpp::RcppArmadillo::sample(idx, size, true, probs);
}

arma::mat mvrnorm(arma::uword n, const arma::vec mu, const arma::mat sigma_chol){
  int ncols = sigma_chol.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * sigma_chol;
}

arma::vec mvrnorm(const arma::vec mu, const arma::mat sigma_chol){
  return mvrnorm(1, mu, sigma_chol).row(0);
}

