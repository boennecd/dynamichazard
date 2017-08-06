#ifndef DMVNRM
#define DMVNRM

#include "../arma_n_rcpp.h"

/*
  Multivariate normal distribution density. See:
    http://gallery.rcpp.org/articles/dmvnorm_arma/

  Assumes that sigma_chol_inv is upper triangular such that:
    sigma = sigma_chol_inv %*% t(sigma_chol_inv)
*/

inline double dmvnrm_log(const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv){
  const double log2pi = std::log(2.0 * M_PI);
  int xdim = x.n_elem;
  double rootisum = arma::sum(log(sigma_chol_inv.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  arma::vec z = sigma_chol_inv.t() * (x - mean);
  return(constants - 0.5 * arma::sum(z % z) + rootisum);
}

#endif
