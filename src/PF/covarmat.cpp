#include "covarmat.h"

#ifdef _OPENMP
class LockGuard {
  omp_lock_t &m_lock;
public:
  explicit LockGuard(omp_lock_t &lock) : m_lock(lock){
    omp_set_lock(&m_lock);
  }

  ~LockGuard() {
    omp_unset_lock(&m_lock);
  }
};
#endif

const arma::mat& covarmat::get_mat(output what) const {
  if(what == e_mat)
    return(*mat_.get());

  bool *this_flag;

  this_flag = is_chol_set.get();
#ifdef _OPENMP
  bool is_computed;
#pragma omp atomic read
  is_computed = *this_flag;
  if(!is_computed){
    LockGuard guard(*lock);
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
    LockGuard guard(*lock);
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
    LockGuard guard(*lock);
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

const arma::mat& covarmat::mat() const  {
  return get_mat(e_mat);
}
const arma::mat& covarmat::chol() const {
  return get_mat(e_chol);
}
const arma::mat& covarmat::chol_inv() const {
  return get_mat(e_chol_inv);
}
const arma::mat& covarmat::inv() const {
  return get_mat(e_inv);
}

covarmat::covarmat(const covarmat &other): covarmat(other.mat()) { }

#ifdef _OPENMP
covarmat::~covarmat(){
  omp_destroy_lock(lock.get());
}
#endif
