#include "sample_funcs.h"
#include <RcppArmadilloExtensions/sample.h>

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
};


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
};

// --------------------------------------- //

arma::mat mvrnorm(arma::uword n, const arma::vec mu, const arma::mat sigma_chol){
  int ncols = sigma_chol.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * sigma_chol;
}

arma::vec mvrnorm(const arma::vec mu, const arma::mat sigma_chol){
  return mvrnorm(1, mu, sigma_chol).row(0);
}

