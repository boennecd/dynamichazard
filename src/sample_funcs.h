#ifndef SAMPLE_FUNCS
#define SAMPLE_FUNCS
#include "arma_n_rcpp.h"

/*
  Samples with replacement using probs. Like running this R-code:
    probs <- runif(10); probs <- probs / sum(probs)
    -1 + sample.int(length(probs), replace = TRUE, prob = probs)

  The second overloads matches with
     n = 100
     probs <- runif(10); probs <- probs / sum(probs)
     -1 + sample.int(length(probs), n, replace = TRUE, prob = probs)

  See http://gallery.rcpp.org/articles/using-the-Rcpp-based-sample-implementation/
*/
arma::uvec sample_indices(arma::vec probs);
arma::uvec sample_indices(int size, arma::vec probs);

/*
  n draws from N(mu, Sigma). See
    http://gallery.rcpp.org/articles/simulate-multivariate-normal/
*/

arma::mat mvrnorm(arma::uword n, const arma::vec mu, const arma::mat sigma_chol);
arma::vec mvrnorm(const arma::vec mu, const arma::mat sigma_chol);

#endif
