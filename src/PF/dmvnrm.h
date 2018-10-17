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

double dmvnrm_log(arma::vec, const arma::mat&);
double dmvnrm_log(const arma::vec&, const arma::vec&, const arma::mat&);

double dmvtrm_log(arma::vec, const arma::mat&, const int);
double dmvtrm_log(
    const arma::vec&, const arma::vec&, const arma::mat&, const int);

#endif
