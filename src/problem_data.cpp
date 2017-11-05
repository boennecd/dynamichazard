#include "problem_data.h"

using ptr_vec = std::unique_ptr<arma::vec>;
using ptr_mat = std::unique_ptr<arma::mat>;

map_res_col ddhazard_data::get_dynamic_coefs(unsigned int t){
  arma::vec temp(a_t_less_s.colptr(t - 1), a_t_less_s.n_rows, false);
  return lp_map(temp);
}

map_res_col ddhazard_data_random_walk::lp_map(arma::vec &a){
  auto span_use = (order == 1) ? arma::span::all : *span_current_cov;
  return map_res_col(a(span_use));
}

map_res_mat ddhazard_data_random_walk::lp_map(arma::mat &M){
  auto span_use = (order == 1) ? arma::span::all : *span_current_cov;
  return map_res_mat(M(span_use, span_use));
}

map_res_col ddhazard_data_random_walk::lp_map_inv(arma::vec &a){
  if(order == 1)
    return map_res_col(a(arma::span::all));

  ptr_vec ptr(new arma::vec(space_dim, arma::fill::zeros));
  arma::vec &out = *ptr.get();
  out(*span_current_cov) = a;

  return map_res_col(out(arma::span::all), ptr);
}

map_res_mat ddhazard_data_random_walk::lp_map_inv(arma::mat &M){
  if(order == 1)
    return map_res_mat(M(arma::span::all, arma::span::all));

  ptr_mat ptr(new arma::mat(space_dim, space_dim, arma::fill::zeros));
  arma::mat &out = *ptr.get();
  out(*span_current_cov, *span_current_cov) = M;

  return map_res_mat(out(arma::span::all, arma::span::all), ptr);
}
