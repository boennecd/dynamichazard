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

#ifdef _OPENMP
#pragma omp parallel
{
#endif

  arma::mat my_ans(n_elem, n_elem, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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
    sym_mat_rank_one_update(weight, state, my_ans);
    sym_mat_rank_one_update(weight, m, my_ans);
    R_BLAS_LAPACK::dger(
      &n_elem /* M */, &n_elem /* N */, &neg_weight /* ALPHA */,
      state.memptr() /* X */, &inc /* INCX*/,
      m.memptr() /* Y */, &inc /* INCY */,
      my_ans.memptr() /* A */, &n_elem /* LDA */);
    R_BLAS_LAPACK::dger(
      &n_elem, &n_elem, &neg_weight,
      m.memptr() /* swapped */, &inc ,
      state.memptr() /* swapped */, &inc,
      my_ans.memptr(), &n_elem);
  }
#ifdef _OPENMP
#pragma omp critical(get_E_x_less_x_less_one_outer_at_one)
{
#endif
  ans += my_ans;
#ifdef _OPENMP
}
}
#endif

  ans += arma::inv(S_inv);

  return ans;
}

static PF_summary_stats_RW compute_summary_stats_first_o_RW(
    const std::vector<cloud> &smoothed_clouds,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0){
  PF_summary_stats_RW ans;
  std::vector<arma::vec> &E_xs = ans.E_xs;
  std::vector<arma::mat> &E_x_less_x_less_one_outers = ans.E_x_less_x_less_one_outers;

  unsigned int n_periods = smoothed_clouds.size();
  unsigned int n_elem = smoothed_clouds[0][0].get_state().n_elem;

  auto it_cl = smoothed_clouds.begin();
  for(unsigned int i = 0; i < n_periods; ++i, ++it_cl){
    const bool is_first = i == 0;
    const bool is_last = i == n_periods - 1;

    unsigned int n_part = it_cl->size();
    arma::vec E_x(n_elem, arma::fill::zeros);
    arma::mat E_x_less_x_less_one_outer(n_elem, n_elem, arma::fill::zeros);

    auto it_p = it_cl->begin();
    for(unsigned int j = 0; j < n_part; ++j, ++it_p){
      double weight = exp(it_p->log_weight);

      E_x += weight * it_p->get_state();
      if(is_last || is_first)
        continue;

      arma::vec inter = it_p->get_state() - it_p->parent->get_state();
      sym_mat_rank_one_update(weight, inter, E_x_less_x_less_one_outer);
    }

    if(is_last){
      --it_cl;
      for(auto it_p = it_cl->begin(); it_p != it_cl->end(); ++it_p){
        arma::vec inter = it_p->child->get_state() - it_p->get_state();
        sym_mat_rank_one_update(exp(it_p->log_weight), inter, E_x_less_x_less_one_outer);
      }

    } else if(is_first){
      E_x_less_x_less_one_outer =
        get_E_x_less_x_less_one_outer_at_one(a_0, Q, Q_0, *it_cl);
    }

    E_xs.push_back(std::move(E_x));
    E_x_less_x_less_one_outer = arma::symmatu(E_x_less_x_less_one_outer);
    E_x_less_x_less_one_outers.push_back(std::move(E_x_less_x_less_one_outer));
  }

  return ans;
}

static PF_summary_stats_RW compute_summary_stats_first_o_RW(
    const std::vector<std::vector<smoother_output::particle_pairs>> &transition_likelihoods,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
    const cloud &first_smoothed_cloud){
  PF_summary_stats_RW ans;
  std::vector<arma::vec> &E_xs = ans.E_xs;
  std::vector<arma::mat> &E_x_less_x_less_one_outers = ans.E_x_less_x_less_one_outers;

  unsigned int n_periods = transition_likelihoods.size();
  unsigned int n_elem = transition_likelihoods[0][0].p->get_state().n_elem;

  auto it_trans = transition_likelihoods.begin();
  for(unsigned int i = 0; i < n_periods; ++i, ++it_trans){
    const bool is_first = i == 0;

    unsigned int n_part = it_trans->size();
    arma::vec E_x(n_elem, arma::fill::zeros);
    arma::mat E_x_less_x_less_one_outer(n_elem, n_elem, arma::fill::zeros);

    auto it_trans_begin = it_trans->begin();
#ifdef _OPENMP
#pragma omp parallel if(!is_first)
{
#endif
  arma::vec my_E_x(n_elem, arma::fill::zeros);
  arma::mat my_E_x_less_x_less_one_outer(n_elem, n_elem, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for(unsigned int j = 0; j < n_part; ++j){
    auto it_elem = it_trans_begin + j;
    const particle *this_p = it_elem->p;
    double weight_outer = exp(this_p->log_weight);
    my_E_x += weight_outer * this_p->get_state();

    if(is_first)
      continue;

    for(auto it_pair = it_elem->transition_pairs.begin();
        it_pair != it_elem->transition_pairs.end(); ++it_pair){
      const particle *pair_p = it_pair->p;
      double weight_inner = exp(this_p->log_weight + it_pair->log_weight);

      arma::vec inter =  this_p->get_state() - pair_p->get_state();
      sym_mat_rank_one_update(weight_inner, inter, my_E_x_less_x_less_one_outer);
    }
  }
#ifdef _OPENMP
#pragma omp critical(compute_summary_stats_first_o_RW_w_tran_like)
{
#endif
  E_x += my_E_x;
  E_x_less_x_less_one_outer += my_E_x_less_x_less_one_outer;
#ifdef _OPENMP
}
}
#endif

if(is_first){
  E_x_less_x_less_one_outer =
    get_E_x_less_x_less_one_outer_at_one(a_0, Q, Q_0, first_smoothed_cloud);
}

E_xs.push_back(std::move(E_x));
E_x_less_x_less_one_outer = arma::symmatu(E_x_less_x_less_one_outer);
E_x_less_x_less_one_outers.push_back(std::move(E_x_less_x_less_one_outer));
  }

  return ans;
}

PF_summary_stats_RW
  compute_summary_stats_first_o_RW
  (const smoother_output &sm_output, const arma::vec &a_0, const arma::mat &Q,
   const arma::mat &Q_0){
    if(sm_output.transition_likelihoods.size() == 0){
      return(compute_summary_stats_first_o_RW(
          sm_output.smoothed_clouds, a_0, Q, Q_0));
    }

    return(compute_summary_stats_first_o_RW(
        sm_output.transition_likelihoods, a_0, Q, Q_0,
        sm_output.smoothed_clouds.front()));
}
