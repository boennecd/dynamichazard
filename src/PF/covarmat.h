#ifndef COVARMAT_H
#define COVARMAT_H
#include <memory>
#include "../arma_n_rcpp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* data class for pre-computed factorization and matrices which are only
 * computed once */
class covarmat {
private:
#ifdef _OPENMP
  std::unique_ptr<omp_lock_t> lock;
#endif

  enum output { e_mat, e_chol, e_chol_inv, e_inv };

  std::unique_ptr<const arma::mat> mat_;
  std::unique_ptr<bool> is_chol_set =
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> chol_;
  std::unique_ptr<bool> is_chol_inv_set=
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> chol_inv_;
  std::unique_ptr<bool> is_inv_set=
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> inv_;

  const arma::mat& get_mat(output) const;

public:
  /* Covar matrix: Q */
  const arma::mat& mat() const;
  /* C in Q = C^\topC in  */
  const arma::mat& chol() const;
  /* C^{-1} */
  const arma::mat& chol_inv() const;
  /* Q^{-1} */
  const arma::mat& inv() const;

  template<typename T>
  covarmat(T Q):
    mat_     (new arma::mat(Q)),
    chol_    (new arma::mat(arma::size(Q), arma::fill::zeros)),
    chol_inv_(new arma::mat(arma::size(Q), arma::fill::zeros)),
    inv_     (new arma::mat(arma::size(Q), arma::fill::zeros)) {
#ifdef _OPENMP
    lock.reset(new omp_lock_t());
    omp_init_lock(lock.get());
#endif
  }

  covarmat(const covarmat&);

#ifdef _OPENMP
  ~covarmat();
#endif
};

#endif
