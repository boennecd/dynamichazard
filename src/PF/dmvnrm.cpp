#include "dmvnrm.h"

double dmvnrm_log(
    arma::vec z /* x - mean */, const arma::mat &sigma_chol_inv){
  const double log2pi = std::log(2.0 * M_PI);
  int xdim = z.n_elem;
  double rootisum = arma::sum(log(sigma_chol_inv.diag()));
  double constants = -((double)xdim/2.0) * log2pi;

  // z = sigma_chol_inv.t() * z;
  int incx = 1;
  char uplo = 'U', trans = 'T', diag = 'N';
  R_BLAS_LAPACK::dtrmv(
    &uplo, &trans, &diag,
    &xdim /* N */, sigma_chol_inv.memptr() /* A */,
    &xdim /* LDA */, z.memptr() /* X */, &incx /* INCX */);

  return(constants - 0.5 * arma::sum(z % z) + rootisum);
}

double dmvnrm_log(
    const arma::vec &x, const arma::vec &mean,
    const arma::mat &sigma_chol_inv){
  return dmvnrm_log(x - mean, sigma_chol_inv);
}



double dmvtrm_log(
    arma::vec z /* x - mean */, const arma::mat &sigma_chol_inv, const int nu){
  int xdim = z.n_elem;
  const double
    dnu = nu, p = xdim,
    rootisum = arma::sum(log(sigma_chol_inv.diag())),
    constants = std::lgamma((p + dnu) * .5) - lgamma(dnu * .5) -
      std::log(dnu * M_PI) * p * .5 + rootisum;

  /* z = sigma_chol_inv.t() * z; */
  int incx = 1;
  char uplo = 'U', trans = 'T', diag = 'N';
  R_BLAS_LAPACK::dtrmv(
    &uplo, &trans, &diag,
    &xdim /* N */, sigma_chol_inv.memptr() /* A */,
    &xdim /* LDA */, z.memptr() /* X */, &incx /* INCX */);

  return constants - (dnu + p) * .5 *  std::log1p(arma::sum(z % z) / dnu);
}

double dmvtrm_log(
    const arma::vec &x, const arma::vec &mean,
    const arma::mat &sigma_chol_inv, const int nu){
  return dmvtrm_log(x - mean, sigma_chol_inv, nu);
}
