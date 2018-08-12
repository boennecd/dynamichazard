#include "est_params.h"
#include "../arma_BLAS_LAPACK.h"
#include "../R_BLAS_LAPACK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static PF_summary_stats compute_PF_summary_stats(
    const smoother_output::trans_like_obj &transition_likelihoods,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
    const cloud &first_smoothed_cloud, const arma::mat F,
    const bool do_use_F, const bool do_compute_E_x){
  PF_summary_stats ans;
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

    arma::vec &E_x = E_xs[i];
    E_x.zeros(n_elem);
    arma::mat &E_x_less_x_less_one_outer =
      E_x_less_x_less_one_outers[i];
    E_x_less_x_less_one_outer.zeros(n_elem, n_elem);

    auto it_trans_begin = it_trans->begin();
    for(auto it_elem = it_trans_begin; it_elem != it_trans->end(); ++it_elem){
      const particle *this_p = it_elem->p;
      if(do_compute_E_x)
        E_x += exp(it_elem->log_weight) * this_p->get_state();

      for(auto it_pair = it_elem->transition_pairs.begin();
          it_pair != it_elem->transition_pairs.end(); ++it_pair){
        const particle *pair_p = it_pair->p;
        double weight_inner = exp(it_elem->log_weight + it_pair->log_weight);

        arma::vec inter = do_use_F ?
          arma::vec(this_p->get_state() - F * pair_p->get_state()) :
          arma::vec(this_p->get_state() -     pair_p->get_state());

        sym_mat_rank_one_update(
          weight_inner, inter, E_x_less_x_less_one_outer);
      }
    }

    E_x_less_x_less_one_outer = arma::symmatu(E_x_less_x_less_one_outer);
  }

  return ans;
}

PF_summary_stats
  compute_PF_summary_stats
  (const smoother_output &sm_output, const arma::vec &a_0, const arma::mat &Q,
   const arma::mat &Q_0, const arma::mat F, const bool do_use_F,
   const bool do_compute_E_x){

    std::shared_ptr<smoother_output::trans_like_obj> trans_obj_ptr =
      sm_output.get_transition_likelihoods(true);

    return(compute_PF_summary_stats(
        *trans_obj_ptr, a_0, Q, Q_0, sm_output.smoothed_clouds.front(),
        F, do_use_F, do_compute_E_x));
}
