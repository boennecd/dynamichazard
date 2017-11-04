#include "problem_data.h"

coefs problem_data_random_walk::get_coefs(unsigned int t){
  arma::vec tmp(a_t_less_s.colptr(t - 1), a_t_less_s.n_rows, false);
  return { get_coefs(tmp), fixed_parems };
}

arma::vec problem_data_random_walk::get_coefs(arma::vec &a){
  return a(*span_current_cov);
}
