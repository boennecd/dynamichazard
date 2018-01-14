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
arma::uvec sample_indices(arma::vec &probs);
arma::uvec sample_indices(const arma::uword size, arma::vec &probs);

/*
  Systematic resampling function from from:
    Kitagawa, G. (1996). Monte Carlo filter and smoother for non-Gaussian nonlinear state space models. Journal of computational and graphical statistics, 5(1), 1-25.
*/
arma::uvec systematic_resampling(arma::vec &probs);
arma::uvec systematic_resampling(const arma::uword size, arma::vec &probs);

/*
  n draws from N(mu, Sigma). See
    http://gallery.rcpp.org/articles/simulate-multivariate-normal/
*/

arma::mat mvrnorm(const int, const arma::mat&);
arma::mat mvrnorm(const int, const arma::vec&, const arma::mat&);
arma::vec mvrnorm(const arma::vec&, const arma::mat&);
arma::vec mvrnorm(const arma::mat&);

#endif
