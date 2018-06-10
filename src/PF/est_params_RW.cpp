#include "PF_utils.h"
#include "../arma_BLAS_LAPACK.h"
#include "../R_BLAS_LAPACK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static arma::mat get_E_x_less_x_less_one_outer_at_one(
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
    const cloud &cl){
  const int n_elem = a_0.n_elem;
  const unsigned int n_particles = cl.size();
  auto cl_begin = cl.begin();

  arma::mat ans(n_elem, n_elem, arma::fill::zeros);
  const arma::mat S_inv = arma::inv(Q) + arma::inv(Q_0);
  const arma::vec a_0_term = solve(Q_0, a_0);

  const arma::mat Q_chol = arma::chol(Q);
  const arma::mat S_inv_chol = arma::chol(S_inv);

  for(unsigned int i = 0; i < n_particles; ++i){
    auto it_p = cl_begin + i;
    const arma::vec &state = it_p->get_state();
    arma::vec m = a_0_term + solve_w_precomputed_chol(Q_chol, state);
    m = solve_w_precomputed_chol(S_inv_chol, m);
    double weight = exp(it_p->log_weight);
    double neg_weight = -weight;
    const static int inc = 1;

    // Could use dsyrk and dsyr2k
    // This parts only sets the upper triangular part of the matrix
    sym_mat_rank_one_update(weight, state, ans);
    sym_mat_rank_one_update(weight, m, ans);
    R_BLAS_LAPACK::dger(
      &n_elem /* M */, &n_elem /* N */, &neg_weight /* ALPHA */,
      state.memptr() /* X */, &inc /* INCX*/,
      m.memptr() /* Y */, &inc /* INCY */,
      ans.memptr() /* A */, &n_elem /* LDA */);
    R_BLAS_LAPACK::dger(
      &n_elem, &n_elem, &neg_weight,
      m.memptr() /* swapped */, &inc ,
      state.memptr() /* swapped */, &inc,
      ans.memptr(), &n_elem);
  }

  ans += arma::inv(S_inv);

  return ans;
}

static PF_summary_stats_RW compute_summary_stats_first_o_RW(
    const smoother_output::trans_like_obj &transition_likelihoods,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
    const cloud &first_smoothed_cloud){
  PF_summary_stats_RW ans;
  unsigned int n_periods = transition_likelihoods.size();
  unsigned int n_elem = transition_likelihoods[0][0].p->get_state().n_elem;
  std::vector<arma::vec> &E_xs = ans.E_xs;
  std::vector<arma::mat> &E_x_less_x_less_one_outers =
    ans.E_x_less_x_less_one_outers;
  E_xs.resize(n_periods);
  E_x_less_x_less_one_outers.resize(n_periods);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned int i = 0; i < n_periods; ++i){
    auto it_trans = transition_likelihoods.begin() + i;
    const bool is_first = i == 0;

    arma::vec &E_x = E_xs[i];
    E_x.zeros(n_elem);
    arma::mat &E_x_less_x_less_one_outer =
      E_x_less_x_less_one_outers[i];
    E_x_less_x_less_one_outer.zeros(n_elem, n_elem);

    auto it_trans_begin = it_trans->begin();
    for(auto it_elem = it_trans_begin; it_elem != it_trans->end(); ++it_elem){
      const particle *this_p = it_elem->p;
      E_x += exp(it_elem->log_weight) * this_p->get_state();

      if(is_first)
        continue;

      for(auto it_pair = it_elem->transition_pairs.begin();
          it_pair != it_elem->transition_pairs.end(); ++it_pair){
        const particle *pair_p = it_pair->p;
        double weight_inner = exp(it_elem->log_weight + it_pair->log_weight);

        arma::vec inter =  this_p->get_state() - pair_p->get_state();
        sym_mat_rank_one_update(
          weight_inner, inter, E_x_less_x_less_one_outer);
      }
    }

    if(is_first){
      E_x_less_x_less_one_outer =
        get_E_x_less_x_less_one_outer_at_one(
          a_0, Q, Q_0, first_smoothed_cloud);
    }

    E_x_less_x_less_one_outer = arma::symmatu(E_x_less_x_less_one_outer);
  }

  return ans;
}

PF_summary_stats_RW
  compute_summary_stats_first_o_RW
  (const smoother_output &sm_output, const arma::vec &a_0, const arma::mat &Q,
   const arma::mat &Q_0){

    std::shared_ptr<smoother_output::trans_like_obj> trans_obj_ptr =
      sm_output.get_transition_likelihoods(true);

    return(compute_summary_stats_first_o_RW(
        *trans_obj_ptr, a_0, Q, Q_0, sm_output.smoothed_clouds.front()));
}
