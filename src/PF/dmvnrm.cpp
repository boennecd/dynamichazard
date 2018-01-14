#include "dmvnrm.h"

double dmvnrm_log(
    arma::vec z /* x - mean */, const arma::mat &sigma_chol_inv){
  const double log2pi = std::log(2.0 * M_PI);
  int xdim = z.n_elem;
  double rootisum = arma::sum(log(sigma_chol_inv.diag()));
  double constants = -((double)xdim/2.0) * log2pi;

  // z = sigma_chol_inv.t() * z;
  int incx = 1;
  R_BLAS_LAPACK::dtrmv(
    "U" /* UPLO */, "T" /* TRANS */, "N" /* diag */,
    &xdim /* N */, sigma_chol_inv.memptr() /* A */,
    &xdim /* LDA */, z.memptr() /* X */, &incx /* INCX */);
    return(constants - 0.5 * arma::sum(z % z) + rootisum);
}

double dmvnrm_log(
    const arma::vec &x, const arma::vec &mean,
    const arma::mat &sigma_chol_inv){
  return dmvnrm_log(x - mean, sigma_chol_inv);
}
