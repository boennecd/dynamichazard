#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"
#include "arma_BLAS_LAPACK.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Util function
static inline unsigned int get_cloud_idx(
    const particle *const p){
  return ((p) ? p->get_cloud_idx() + 1  /* want non-zero based */ : 0);
}

static inline const particle* get_cloud_idx(
    const cloud &cl, const unsigned int r_idx){
  return r_idx > 0 ? &cl[r_idx - 1]  /* not zero based */ : nullptr;
}

/* Function to turn a clouds of particles into an Rcpp::List */
Rcpp::List get_rcpp_list_from_cloud(
    const std::vector<cloud> &clouds, const bool reverse,
    unsigned int state_dim, const PF_data *data){
  if(data && data->debug > 0)
    data->log(1) << "Creating output list";

  auto n_clouds = clouds.size();
  Rcpp::List ans(n_clouds);
  std::vector<cloud>::const_iterator it =
    reverse ? clouds.end() - 1 : clouds.begin();
  for(std::vector<cloud>::size_type j = 0;
      j < n_clouds;
      ++j, it += 1 - 2 * reverse){
    auto n_elem = it->size();
    arma::uvec parent_idx(n_elem);
    arma::uvec child_idx(n_elem);
    arma::vec log_unnormalized_weights(n_elem);
    arma::vec weights(n_elem);
    arma::mat states(state_dim, n_elem);

    auto idx_pr = parent_idx.begin();
    auto idx_ch = child_idx.begin();
    auto un_w = log_unnormalized_weights.begin();
    auto w = weights.begin();
    auto pr = it->begin();
    for(arma::uword i = 0;
        i < n_elem;
        ++i, ++w, ++idx_pr, ++idx_ch, ++pr, ++un_w){
      *idx_pr = get_cloud_idx(pr->parent);
      *idx_ch = get_cloud_idx(pr->child);

      *w = exp(pr->log_weight);
      *un_w = pr->log_unnormalized_weight;

      states.col(i) = pr->get_state();
    }

    ans[j] =
      Rcpp::List::create(
        Rcpp::Named("parent_idx") = Rcpp::wrap(parent_idx),
        Rcpp::Named("child_idx") = Rcpp::wrap(child_idx),
        Rcpp::Named("weights") = Rcpp::wrap(weights),
        Rcpp::Named("states") = Rcpp::wrap(states),
        Rcpp::Named("log_unnormalized_weights") = Rcpp::wrap(log_unnormalized_weights)
      );
  }

  return ans;
}

Rcpp::List get_rcpp_list_from_cloud(
    const smoother_output &sm_output, const PF_data *data){
  unsigned int state_dim = sm_output.forward_clouds[0][0].get_state().n_elem;

  // Set-up the transition_likelihoods element
  Rcpp::List transitions; // List with final output
  auto &transition_likelihoods = sm_output.transition_likelihoods;
  if(transition_likelihoods.size() > 0){
    transitions = Rcpp::List(transition_likelihoods.size());

    auto it_trans = transitions.begin();
    // Time loop
    for(auto it = transition_likelihoods.begin();
        it != transition_likelihoods.end(); ++it, ++it_trans){
      Rcpp::List ans_inner(it->size());

      // Loop over first index of pair of particles
      auto it_ans = ans_inner.begin();
      for(auto it_inner = it->begin();
          it_inner != it->end(); ++it_inner, ++it_ans){
        const unsigned int n_inner = it_inner->transition_pairs.size();
        arma::uvec fw_idx(n_inner);
        arma::vec weights(n_inner);

        auto it_fw_idx = fw_idx.begin();
        auto it_w = weights.begin();
        auto it_data = it_inner->transition_pairs.begin();
        // Loop over second index of pair of particles
        for(unsigned int j = 0; j < n_inner; ++j, ++it_fw_idx, ++it_w, ++it_data){
          *it_fw_idx = get_cloud_idx(it_data->p);
          *it_w = exp(it_data->log_weight);
        }

        *it_ans = Rcpp::List::create(
          Rcpp::Named("fw_idx") = Rcpp::wrap(std::move(fw_idx)),
          Rcpp::Named("weights") = Rcpp::wrap(std::move(weights)),
          Rcpp::Named("bw_idx") = get_cloud_idx(it_inner->p)
        );
      }

      *it_trans = std::move(ans_inner);
    }
  }

  // Create output list
  Rcpp::List ans = Rcpp::List::create(
    Rcpp::Named("forward_clouds") = get_rcpp_list_from_cloud(
      sm_output.forward_clouds, false, state_dim, data),
    Rcpp::Named("backward_clouds") = get_rcpp_list_from_cloud(
      sm_output.backward_clouds, true, state_dim, data),
    Rcpp::Named("smoothed_clouds") = get_rcpp_list_from_cloud(
      sm_output.smoothed_clouds, false, state_dim, data),
    Rcpp::Named("transition_likelihoods") = std::move(transitions));

  ans.attr("class") = "PF_clouds";

  return ans;
}

/* Function to turn an Rcpp::List into a clouds of particles */
template<bool is_smooth, const bool reverse>
static std::vector<cloud> get_clouds_from_rcpp_list_util
  (const Rcpp::List &rcpp_list,
   const std::vector<cloud> *fw = nullptr,
   const std::vector<cloud> *bw = nullptr)
  {
    unsigned int n_periods = rcpp_list.size();
    std::vector<cloud> ans(n_periods);

    auto it_ans = ans.begin();
    auto it_ans_prev = it_ans - 1;
    for(unsigned int i = reverse ? n_periods - 1 : 0;
        reverse ? true /* we break later */ : i < n_periods;
        i += 1 - 2 * reverse, ++it_ans, ++it_ans_prev){
      Rcpp::List cloud_list(rcpp_list[i]);

      arma::uvec parent_idx = Rcpp::as<arma::uvec>(cloud_list["parent_idx"]);
      arma::uvec child_idx = Rcpp::as<arma::uvec>(cloud_list["child_idx"]);
      arma::vec weights = Rcpp::as<arma::vec>(cloud_list["weights"]);
      arma::vec log_unnormalized_weights = Rcpp::as<arma::vec>(
        cloud_list["log_unnormalized_weights"]);
      arma::mat states = Rcpp::as<arma::mat>(cloud_list["states"]);

      auto n_states = weights.n_elem;
      auto it_par = parent_idx.begin();
      auto it_child = child_idx.begin();
      auto it_w = weights.begin();
      auto it_un_w = log_unnormalized_weights.begin();
      for(unsigned j = 0; j < n_states; ++j, ++it_w, ++it_par, ++it_child, ++it_un_w){
        const particle *parent;
        if (*it_par == 0){
          parent = nullptr;

        } else if(is_smooth){
          parent = get_cloud_idx((*fw)[i /* starts at time 0 */], *it_par);

        } else
          parent = get_cloud_idx((*it_ans_prev), *it_par);

        const particle *child  =
          /* starts at time 1 and has on more index */
          is_smooth ?
            get_cloud_idx((*bw)[n_periods - (i + 1)], *it_child) : nullptr;

        it_ans->new_particle(states.col(j), parent, child);
        particle &p = it_ans->back();
        p.log_weight = log(*it_w);
        p.log_unnormalized_weight = *it_un_w;
      }

      if(reverse && i == 0) break;
    }

    return ans;
}

smoother_output get_clouds_from_rcpp_list(const Rcpp::List &rcpp_list){
  smoother_output ans;

  // Find clouds
  std::vector<cloud> &fw_clouds = ans.forward_clouds;
  std::vector<cloud> &bw_clouds = ans.backward_clouds;
  std::vector<cloud> &smooth_clouds = ans.smoothed_clouds;

  fw_clouds = get_clouds_from_rcpp_list_util<false, false>(rcpp_list["forward_clouds"]);
  bw_clouds = get_clouds_from_rcpp_list_util<false, true>(rcpp_list["backward_clouds"]);
  smooth_clouds = get_clouds_from_rcpp_list_util<true, false>(
    rcpp_list["smoothed_clouds"], &fw_clouds, &bw_clouds);

  // Set transition likelihoods
  std::vector<std::vector<smoother_output::particle_pairs>>
    &transition_likelihoods = ans.transition_likelihoods;
  const Rcpp::List transitions((SEXP)rcpp_list["transition_likelihoods"]);
  if(transitions.size() > 0){
    unsigned int n_periods = transitions.size();
    transition_likelihoods.resize(n_periods);

    auto it_trans = transitions.begin();
    auto it_ans = transition_likelihoods.begin();
    // Time loop
    for(unsigned int i = 0; i < n_periods; ++i, ++it_trans, ++it_ans){
      Rcpp::List input_at_t(*it_trans);
      unsigned int n_elem = input_at_t.size();
      std::vector<smoother_output::particle_pairs> &new_elems = *it_ans;
      new_elems.resize(n_elem);

      cloud &fw_cloud = fw_clouds[i]; // starts at time 0
      cloud &sm_cloud = smooth_clouds[i]; // starts at time 1

      auto it_input_at_t = input_at_t.begin();
      auto it_elems = new_elems.begin();
      // Loop over first index of pair of particles
      for(unsigned int j = 0; j < n_elem; ++j, ++it_input_at_t, ++it_elems){
        smoother_output::particle_pairs &elem = *it_elems;
        std::vector<smoother_output::pair> &transition_pairs =
          elem.transition_pairs;

        Rcpp::List elem_info(*it_input_at_t);
        unsigned int idx = Rcpp::as<unsigned int>(elem_info["bw_idx"]);
        elem.p = get_cloud_idx(sm_cloud, idx);

        arma::uvec fw_idx = Rcpp::as<arma::uvec>(elem_info["fw_idx"]);
        arma::vec weights = Rcpp::as<arma::vec>(elem_info["weights"]);

        unsigned int n_elem_inner = fw_idx.n_elem;
        transition_pairs.resize(n_elem_inner);
        auto it_idx = fw_idx.begin();
        auto it_w = weights.begin();
        auto it_trans_pair = transition_pairs.begin();
        // Loop over second index of pair of particles
        for(unsigned int k = 0; k < n_elem_inner; ++k, ++it_idx, ++it_w, ++it_trans_pair){
          it_trans_pair->p = get_cloud_idx(fw_cloud, *it_idx);
          it_trans_pair->log_weight = log(*it_w);
        }
      }
    }
  }

  return ans;
}

/* --------------------------------------- */

template<
  template <
    template <typename, bool> class,
    template <typename, bool> class,
    class>
  class smoother>
class PF_smooth_T {
public:
  using bootstrap_filter_sm =
    smoother<
      None_AUX_resampler,
      importance_dens_no_y_dependence,
      binary>;

  using PF_w_normal_approx_sm =
    smoother<
      None_AUX_resampler,
      importance_dens_normal_approx_w_cloud_mean,
      binary>;

  using AUX_w_normal_approx_sm =
    smoother<
      AUX_resampler_normal_approx_w_cloud_mean,
      importance_dens_normal_approx_w_cloud_mean,
      binary>;

  using PF_w_particles_sm =
    smoother<
      None_AUX_resampler,
      importance_dens_normal_approx_w_particles,
      binary>;

  using AUX_w_particles_sm =
    smoother<
      AUX_resampler_normal_approx_w_particles,
      importance_dens_normal_approx_w_particles,
      binary>;

  static Rcpp::List compute(
      const PF_data &data, const std::string method){

    /* Get the smoothed particles at time 1, 2, ..., d */
    smoother_output result;
    if (method == "bootstrap_filter") {
      result = bootstrap_filter_sm::compute(data);

    } else if (method == "PF_normal_approx_w_cloud_mean"){
      result = PF_w_normal_approx_sm::compute(data);

    } else if (method == "AUX_normal_approx_w_cloud_mean"){
      result = AUX_w_normal_approx_sm::compute(data);

    } else if (method == "PF_normal_approx_w_particles"){
      result = PF_w_particles_sm::compute(data);

    }  else if (method == "AUX_normal_approx_w_particles"){
      result = AUX_w_particles_sm::compute(data);

    } else {
      std::stringstream stream;
      stream << "method '" << method << "' is not implemented";
      Rcpp::stop(stream.str());
    }

    /* Create output list */
    return(get_rcpp_list_from_cloud(result, &data));
  }
};

// [[Rcpp::export]]
Rcpp::List PF_smooth(
    const int n_fixed_terms_in_state_vec,
    arma::mat &X,
    arma::mat &fixed_terms,
    const arma::vec &tstart,
    const arma::vec &tstop,
    const arma::colvec &a_0,
    arma::mat &Q_0,
    arma::mat &Q,
    const arma::mat Q_tilde,
    const Rcpp::List &risk_obj,
    const arma::mat &F,
    const int n_max,
    const int order,
    const int n_threads,
    const int N_fw_n_bw,
    const int N_smooth,
    Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
    const int debug,
    const int N_first,
    const std::string method,
    const std::string smoother){
  const arma::ivec is_event_in_bin = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  PF_data data(
      n_fixed_terms_in_state_vec,
      X,
      fixed_terms,
      tstart,
      tstop,
      is_event_in_bin,
      a_0,
      Q_0,
      Q,
      risk_obj,
      F,
      n_max,
      order,
      n_threads,
      Q_tilde,
      N_fw_n_bw,
      N_smooth,
      forward_backward_ESS_threshold,
      debug,
      N_first);

  Rcpp::List ans;

  if(smoother == "Fearnhead_O_N"){
    ans =
      PF_smooth_T<PF_smoother_Fearnhead_O_N>::compute(data, method);

  } else if (smoother == "Brier_O_N_square"){
    ans =
      PF_smooth_T<PF_smoother_Brier_O_N_square>::compute(data, method);

  } else {
    std::stringstream stream;
    stream << "smoother '" << smoother << "' is not implemented";
    Rcpp::stop(stream.str());

  }

  return(ans);
}

/* --------------------------------------- */

struct PF_summary_stats {
  /* E[X_t] */
  std::vector<arma::vec> E_xs;
  /* E[(X_t - FX_{t-1})(X_t - FX_{t-1})^T] */
  std::vector<arma::mat> E_x_less_x_less_one_outers;
};

static arma::mat get_E_x_less_x_less_one_outer_at_one(
  const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
  const cloud &cl){
  // TODO: make chol and reduce computational cost
  const arma::uword n_elem = a_0.n_elem;
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
    m *= sqrt(weight);
    arma::vec tmp = sqrt(weight) * state;

    // TODO: use rank one updates
    my_ans += tmp * tmp.t()
      - m * tmp.t() - tmp * m.t()
      + m * m.t();
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

static PF_summary_stats compute_summary_stats(
    const std::vector<cloud> &smoothed_clouds,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0){
  PF_summary_stats ans;
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

static PF_summary_stats compute_summary_stats(
    const std::vector<std::vector<smoother_output::particle_pairs>> &transition_likelihoods,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
    const cloud &first_smoothed_cloud){
  PF_summary_stats ans;
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
#pragma omp parallel
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

      if(is_first /* TODO: what to do at time 1 */)
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
#pragma omp critical(compute_summary_stats_w_tran_like)
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

static PF_summary_stats compute_summary_stats(
    const smoother_output sm_output, unsigned int n_threads,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0){
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  if(sm_output.transition_likelihoods.size() == 0){
    return(compute_summary_stats(
        sm_output.smoothed_clouds, a_0, Q, Q_0));
  }

  return(compute_summary_stats(
      sm_output.transition_likelihoods, a_0, Q, Q_0,
      sm_output.smoothed_clouds.front()));
}

// [[Rcpp::export]]
Rcpp::List compute_summary_stats(
    const Rcpp::List &rcpp_list, unsigned int n_threads,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0){
  PF_summary_stats stats;
  {
    auto sm_output = get_clouds_from_rcpp_list(rcpp_list);
    stats = compute_summary_stats(sm_output, n_threads, a_0, Q, Q_0);
  }

  unsigned int n_periods = stats.E_xs.size();
  Rcpp::List ans(n_periods);
  for(unsigned int i = 0; i < n_periods; ++i){
    ans[i] = Rcpp::List::create(
      Rcpp::Named("E_xs") = Rcpp::wrap(stats.E_xs[i]),
      Rcpp::Named("E_x_less_x_less_one_outers") =
        Rcpp::wrap(stats.E_x_less_x_less_one_outers[i])
    );
  }

  return ans;
}


