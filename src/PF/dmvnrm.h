#ifndef DMVNRM
#define DMVNRM

#include "../arma_n_rcpp.h"
#include "../R_BLAS_LAPACK.h"

/*
  Multivariate normal distribution density. See:
    http://gallery.rcpp.org/articles/dmvnorm_arma/

  Assumes that sigma_chol_inv is upper triangular such that:
    sigma^(-1) = sigma_chol_inv %*% t(sigma_chol_inv)
*/

inline double dmvnrm_log(const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv){
  const double log2pi = std::log(2.0 * M_PI);
  int xdim = x.n_elem;
  double rootisum = arma::sum(log(sigma_chol_inv.diag()));
  double constants = -((double)xdim/2.0) * log2pi;

  // z = sigma_chol_inv.t() * (x - mean);
  arma::vec z = x - mean;
  int incx = 1;
  R_BLAS_LAPACK::dtrmv(
    "U" /* UPLO */, "T" /* TRANS */, "N" /* diag */,
    &xdim /* N */, sigma_chol_inv.memptr() /* A */,
    &xdim /* LDA */, z.memptr() /* X */, &incx /* INCX */);
  return(constants - 0.5 * arma::sum(z % z) + rootisum);
}

#endif
