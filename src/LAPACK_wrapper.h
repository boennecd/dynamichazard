#ifndef LAPACK_WRAPPER
#define LAPACK_WRAPPER
#include "arma_n_rcpp.h"

void symmetric_mat_chol(const arma::mat&, arma::mat &);

void chol_rank_one_update(arma::mat&, arma::vec);

void square_tri_inv(const arma::mat&, arma::mat&);

//void square_tri_vec_prod(arma::mat&, const arma::vec&, arma::vec&);

#endif
