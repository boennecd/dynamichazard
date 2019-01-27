#ifndef COVARMAT_H
#define COVARMAT_H

#ifdef _OPENMP
#include <omp.h>
#endif

/* data class for pre-computed factorization and matrices which are only
 * computed once */
class covarmat {
private:
#ifdef _OPENMP
  std::unique_ptr<omp_lock_t> lock =
    std::unique_ptr<omp_lock_t>((new omp_lock_t()));
  class LockGuard {
  public:
    explicit LockGuard(omp_lock_t& lock) : m_lock(lock){
      omp_set_lock(&m_lock);
    }

    ~LockGuard() {
      omp_unset_lock(&m_lock);
    }

  private:
    omp_lock_t& m_lock;
  };
#endif

  enum output { e_mat, e_chol, e_chol_inv, e_inv };

  std::unique_ptr<const arma::mat> mat_;
  std::unique_ptr<bool> is_chol_set =
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> chol_;
  std::unique_ptr<bool> is_chol_inv_set =
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> chol_inv_;
  std::unique_ptr<bool> is_inv_set =
    std::unique_ptr<bool>(new bool(false));
  std::unique_ptr<arma::mat> inv_;

  const arma::mat& get_mat(output what) const {
    if(what == e_mat)
      return(*mat_.get());

    bool is_computed, *this_flag;

    this_flag = is_chol_set.get();
#ifdef _OPENMP
#pragma omp atomic read
    is_computed = *this_flag;
    if(!is_computed){
      LockGuard guard(*lock.get());
#pragma omp atomic read
      is_computed = *this_flag;
      if(!is_computed){
        *chol_.get() += arma::chol(*mat_.get());
#pragma omp atomic write
        *this_flag = true;
      }
    }
#else
    if(!*this_flag){
      *chol_.get() += arma::chol(*mat_.get());
      *this_flag = true;
    }
#endif
    if(what == e_chol)
      return(*chol_.get());


    this_flag = is_chol_inv_set.get();
#ifdef _OPENMP
#pragma omp atomic read
    is_computed = *this_flag;
    if(!is_computed){
      LockGuard guard(*lock.get());
#pragma omp atomic read
      is_computed = *this_flag;
      if(!is_computed){
        *chol_inv_.get() += arma::inv(arma::trimatu(*chol_.get()));
#pragma omp atomic write
        *this_flag = true;
      }
    }
#else
    if(!*this_flag){
      *chol_inv_.get() += arma::inv(arma::trimatu(*chol_.get()));
      *this_flag = true;
    }
#endif
    if(what == e_chol_inv)
      return(*chol_inv_.get());


    this_flag = is_inv_set.get();
#ifdef _OPENMP
#pragma omp atomic read
    is_computed = *this_flag;
    if(!is_computed){
      LockGuard guard(*lock.get());
#pragma omp atomic read
      is_computed = *this_flag;
      if(!is_computed){
        *inv_.get() += mat_->i();
#pragma omp atomic write
        *this_flag = true;
      }
    }
#else
    if(!*this_flag){
      *inv_.get() += mat_->i();
      *this_flag = true;
    }
#endif
    return(*inv_.get());
  }

public:
  /* Covar matrix: Q */
  const arma::mat& mat() const  {
    return get_mat(e_mat);
  }
  /* C in Q = C^\topC in  */
  const arma::mat& chol() const {
    return get_mat(e_chol);
  }
  /* C^{-1} */
  const arma::mat& chol_inv() const {
    return get_mat(e_chol_inv);
  }
  /* Q^{-1} */
  const arma::mat& inv() const {
    return get_mat(e_inv);
  }

  template<typename T>
  covarmat(T Q):
    mat_     (new arma::mat(Q)),
    chol_    (new arma::mat(arma::size(Q), arma::fill::zeros)),
    chol_inv_(new arma::mat(arma::size(Q), arma::fill::zeros)),
    inv_     (new arma::mat(arma::size(Q), arma::fill::zeros)) {
#ifdef _OPENMP
    omp_init_lock(lock.get());
#endif
  }

  covarmat(const covarmat &other): covarmat(other.mat()) { }

  ~covarmat(){
#ifdef _OPENMP
    omp_destroy_lock(lock.get());
#endif
  }
};

#endif
