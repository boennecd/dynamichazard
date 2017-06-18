#include "arma_n_rcpp.h"

void symmetric_mat_chol(const arma::mat&, arma::mat &);

void chol_rank_one_update(arma::mat&, arma::vec);

void square_tri_inv(const arma::mat&, arma::mat&);

void tri_mat_times_vec(arma::mat&, const arma::vec&, arma::vec&, bool);

void tri_mat_times_vec(arma::mat&, arma::vec&, bool);

void sym_mat_rank_one_update(const double, const arma::vec&, arma::mat&);

arma::vec sym_mat_times_vec(const arma::vec&, const arma::mat&);
