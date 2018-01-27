#ifndef ARMA_BLAS_LAPACK
#define ARMA_BLAS_LAPACK

#include "arma_n_rcpp.h"
#include <memory>

void symmetric_mat_chol(const arma::mat&, arma::mat &);

void chol_rank_one_update(arma::mat&, arma::vec);

void square_tri_inv(const arma::mat&, arma::mat&);

void tri_mat_times_vec(const arma::mat&, const arma::vec&, arma::vec&, bool);

void tri_mat_times_vec(const arma::mat&, arma::vec&, bool);

void sym_mat_rank_one_update(const double, const arma::vec&, arma::mat&);

arma::vec sym_mat_times_vec(const arma::vec&, const arma::mat&);

arma::mat out_mat_prod(const arma::mat &);

template<class T>
arma::vec solve_w_precomputed_chol(const arma::mat&, const T&);

class LU_factorization {
  const int M;
  const int N;
  const bool has_elem;
  std::unique_ptr<double []> A;
  std::unique_ptr<int []> IPIV;

public:
  LU_factorization(const arma::mat&);
  arma::mat solve(const arma::mat&, const bool transpose = false) const;
  arma::vec solve(const arma::vec&, const bool transpose = false) const;
};

class QR_factorization {
  const int M;
  const int N;
  std::unique_ptr<double []> qr;
  int rank;
  std::unique_ptr<double []> qraux;
  std::unique_ptr<int []> pivot_;

public:
  QR_factorization(const arma::mat&);
  arma::mat qy(const arma::mat&, const bool transpose = false) const;
  arma::vec qy(const arma::vec&, const bool transpose = false) const;
  arma::mat R() const;
  arma::uvec pivot() const;
};

#endif
