#include "sample_funcs.h"
#include <RcppArmadilloExtensions/sample.h>
#include "R_BLAS_LAPACK.h"

arma::uvec sample_indices(arma::vec &probs){
  return sample_indices(probs.n_elem, probs);
}

arma::uvec sample_indices(const arma::uword size, arma::vec &probs){
  arma::uvec idx = arma::linspace<arma::uvec>(0, probs.n_elem - 1, probs.n_elem);

  return Rcpp::RcppArmadillo::sample(idx, size, true, probs);
}

// --------------------------------------- //

arma::uvec systematic_resampling(arma::vec &probs){
  return(systematic_resampling(probs.n_elem, probs));
}


arma::uvec systematic_resampling(const arma::uword size, arma::vec &probs){
  arma::uvec ans(size);

  double U_delta = 1. / size;
  double U_i = Rcpp::as<double>(Rcpp::runif(1, 0, U_delta));

  auto a = ans.begin();
  auto pr = probs.begin();
  double sum = *(pr++);
  arma::uword pr_idx = 0;
  for(arma::uword i = 0; i < size; ++i, ++a, U_i += U_delta){
    while(true){
      if(pr == probs.end())
        break;

      if(sum <= U_i){
        sum += *(pr++);
        ++pr_idx;
        continue;
      }

      break;
    }

    *a = pr_idx;
  }

  return ans;
}

// --------------------------------------- //

arma::mat mvrnorm(const int m, const arma::mat &sigma_chol){
  const int n = sigma_chol.n_cols;

  Rcpp::NumericVector rng_draw = Rcpp::rnorm(m * n);
  arma::mat Y(&rng_draw[0], m, n, false /* don't copy */);

  // Y <-- Y * chol(Sigma)
  const double alpha = 1.;
  char side = 'R', uplo = 'U', transa = 'N', diag = 'N';
  R_BLAS_LAPACK::dtrmm(
    &side, &uplo, &transa, &diag,
    &m /* M */, &n /* N */, &alpha /* ALPHA */,
    sigma_chol.memptr() /* A */, &n /* LDA */,
    Y.memptr() /* B */, &m /* LDB */);

  return Y.t();
}

arma::mat mvrnorm(
    const int m, const arma::vec &mu, const arma::mat &sigma_chol){
  return arma::repmat(mu, 1, m) + mvrnorm(m, sigma_chol);
}

arma::vec mvrnorm(const arma::vec &mu, const arma::mat &sigma_chol){
  return mvrnorm(1, mu, sigma_chol).col(0);
}

arma::vec mvrnorm(const arma::mat &sigma_chol){
  return mvrnorm(1, sigma_chol).col(0);
}

// --------------------------------------- //

arma::mat mvtrnorm(const int m, const arma::mat &sigma_chol, const int nu){
  arma::mat Y = mvrnorm(m, sigma_chol);
  arma::vec chi_draw = Rcpp::as<arma::vec>(Rcpp::rchisq(m, nu));

  double *u = chi_draw.begin();
  for(int i = 0; i < m; ++i, ++u)
    /* see https://stats.stackexchange.com/a/68490/81865 and
     * https://en.wikipedia.org/wiki/Multivariate_t-distribution#Definition */
    Y.col(i) /= std::sqrt(*u / nu);

  return Y;
}

arma::mat mvtrnorm(
    const int m, const arma::vec &mu, const arma::mat &sigma_chol,
    const int nu)
  {
    return arma::repmat(mu, 1, m) + mvtrnorm(m, sigma_chol, nu);
  }

arma::vec mvtrnorm(
    const arma::vec &mu, const arma::mat &sigma_chol, const int nu)
  {
    return mvtrnorm(1, mu, sigma_chol, nu).col(0);
  }

arma::vec mvtrnorm(const arma::mat &sigma_chol, const int nu){
  return mvtrnorm(1, sigma_chol, nu).col(0);
}
