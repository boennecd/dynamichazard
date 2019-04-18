#include "PF_utils.h"

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
    const unsigned int state_dim, const PF_data *data){
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
    arma::vec log_likelihood_term(n_elem);
    arma::vec weights(n_elem);
    arma::mat states(state_dim, n_elem);

    auto idx_pr = parent_idx.begin();
    auto idx_ch = child_idx.begin();
    auto un_w = log_likelihood_term.begin();
    auto w = weights.begin();
    auto pr = it->begin();
    for(arma::uword i = 0;
        i < n_elem;
        ++i, ++w, ++idx_pr, ++idx_ch, ++pr, ++un_w){
      *idx_pr = get_cloud_idx(pr->parent);
      *idx_ch = get_cloud_idx(pr->child);

      *w = exp(pr->log_weight);
      *un_w = pr->log_likelihood_term;

      states.col(i) = pr->get_state();
    }

    ans[j] =
      Rcpp::List::create(
        Rcpp::Named("parent_idx") = Rcpp::wrap(parent_idx),
        Rcpp::Named("child_idx") = Rcpp::wrap(child_idx),
        Rcpp::Named("weights") = Rcpp::wrap(weights),
        Rcpp::Named("states") = Rcpp::wrap(states),
        Rcpp::Named("log_likelihood_term") =
          Rcpp::wrap(log_likelihood_term));
  }

  return ans;
}

Rcpp::List get_rcpp_list_from_cloud(
    const smoother_output &sm_output, const PF_data *data){
  unsigned int state_dim = sm_output.forward_clouds[0][0].get_state().n_elem;

  // Set-up the transition_likelihoods element
  Rcpp::List transitions; // List with final output
  std::shared_ptr<smoother_output::trans_like_obj> trans_ptr =
    sm_output.get_transition_likelihoods(false);
  smoother_output::trans_like_obj &transition_likelihoods = *trans_ptr;
  if(transition_likelihoods.size() > 0){
    transitions = Rcpp::List(transition_likelihoods.size());

    auto it_trans = transitions.begin();
    // Time loop
    for(auto it = transition_likelihoods.begin();
        it != transition_likelihoods.end(); ++it, ++it_trans){
      unsigned int n_bw = it->size();
      // assume equal number of pair for each element
      unsigned int n_fw = it->front().transition_pairs.size();
      arma::uvec bw_idx(n_bw);
      arma::umat fw_idx(n_fw, n_bw);
      arma::mat weights(n_fw, n_bw);

      // Loop over first index of pair of particles
      auto it_w_begin = weights.begin();
      auto it_fw_begin = fw_idx.begin();
      auto it_bw_begin = bw_idx.begin();
      auto it_inner_begin = it->begin();
#ifdef _OPENMP
      int n_threads_use = std::min(omp_get_max_threads(), (int)n_bw / 200 + 1);
#pragma omp parallel for num_threads(n_threads_use) if(n_threads_use > 1)
#endif
      for(unsigned int j = 0; j < n_bw; ++j){
        auto it_inner = it_inner_begin + j;
        auto it_bw = it_bw_begin + j;
        *it_bw = get_cloud_idx(it_inner->p);

        auto it_data = it_inner->transition_pairs.begin();
        auto it_w = it_w_begin + j * n_fw;
        auto it_fw = it_fw_begin + j * n_fw;
        // Loop over second index of pair of particles
        for(unsigned int i = 0; i < n_fw; ++i, ++it_fw, ++it_w, ++it_data){
          *it_fw = get_cloud_idx(it_data->p);
          *it_w = exp(it_data->log_weight);
        }
      }

      *it_trans = Rcpp::List::create(
        Rcpp::Named("fw_idx") = Rcpp::wrap(std::move(fw_idx)),
        Rcpp::Named("weights") = Rcpp::wrap(std::move(weights)),
        Rcpp::Named("bw_idx") = Rcpp::wrap(std::move(bw_idx))
      );
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
std::vector<cloud> get_cloud_from_rcpp_list
  (const Rcpp::List &rcpp_list,
   const std::vector<cloud> *fw, const std::vector<cloud> *bw)
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
    arma::vec log_likelihood_term = Rcpp::as<arma::vec>(
      cloud_list["log_likelihood_term"]);
    arma::mat states = Rcpp::as<arma::mat>(cloud_list["states"]);

    auto n_states = weights.n_elem;
    auto it_par = parent_idx.begin();
    auto it_child = child_idx.begin();
    auto it_w = weights.begin();
    auto it_un_w = log_likelihood_term.begin();
    for(unsigned j = 0; j < n_states;
        ++j, ++it_w, ++it_par, ++it_child, ++it_un_w){
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
      p.log_likelihood_term = *it_un_w;
    }

    if(reverse && i == 0) break;
  }

  return ans;
}


template
std::vector<cloud> get_cloud_from_rcpp_list<false, false>
(const Rcpp::List&, const std::vector<cloud>*, const std::vector<cloud>*);
template
std::vector<cloud> get_cloud_from_rcpp_list<true, false>
(const Rcpp::List&, const std::vector<cloud>*, const std::vector<cloud>*);
template
  std::vector<cloud> get_cloud_from_rcpp_list<false, true>
  (const Rcpp::List&, const std::vector<cloud>*, const std::vector<cloud>*);

smoother_output get_clouds_from_rcpp_list(const Rcpp::List &rcpp_list){
  smoother_output ans;

  // Find clouds
  std::vector<cloud> &fw_clouds = ans.forward_clouds;
  std::vector<cloud> &bw_clouds = ans.backward_clouds;
  std::vector<cloud> &smooth_clouds = ans.smoothed_clouds;

  fw_clouds =
    get_cloud_from_rcpp_list<false, false>(rcpp_list["forward_clouds"]);
  bw_clouds =
    get_cloud_from_rcpp_list<false, true>(rcpp_list["backward_clouds"]);
  smooth_clouds = get_cloud_from_rcpp_list<true, false>(
    rcpp_list["smoothed_clouds"], &fw_clouds, &bw_clouds);

  // Set transition likelihoods
  std::shared_ptr<smoother_output::trans_like_obj> trans_ptr =
    ans.get_transition_likelihoods(false);
  smoother_output::trans_like_obj &transition_likelihoods = *trans_ptr;
  const Rcpp::List transitions((SEXP)rcpp_list["transition_likelihoods"]);
  if(transitions.size() > 0){
    unsigned int n_periods = transitions.size();
    transition_likelihoods.resize(n_periods);

    auto it_trans = transitions.begin();
    auto it_ans = transition_likelihoods.begin();
    // Time loop
    for(unsigned int i = 0; i < n_periods; ++i, ++it_trans, ++it_ans){
      cloud &fw_cloud = fw_clouds[i]; // starts at time 0
      cloud &sm_cloud = smooth_clouds[i]; // starts at time 1

      const Rcpp::List inner_list(*it_trans);
      const Rcpp::IntegerVector bw_idx((SEXP)inner_list["bw_idx"]);
      const Rcpp::IntegerMatrix fw_idx((SEXP)inner_list["fw_idx"]);
      const Rcpp::NumericMatrix weights((SEXP)inner_list["weights"]);

      unsigned int n_bw = bw_idx.size();
      unsigned int n_fw = fw_idx.nrow();

      std::vector<smoother_output::particle_pairs> &new_elems = *it_ans;
      new_elems.resize(n_bw);

      // Loop over first index of pair of particles
      auto it_elems_begin = new_elems.begin();
      const int *it_bw_begin = &bw_idx[0];
      const int *it_fw_begin = &fw_idx(0,0);
      const double *it_w_begin = &weights(0,0);
#ifdef _OPENMP
      int n_threads_use = std::min((long int)omp_get_max_threads(),
                                   ((int)n_bw - 1L) / 200L + 1L);
#pragma omp parallel for schedule(static) num_threads(n_threads_use) if(n_threads_use > 1)
#endif
      for(unsigned int j = 0; j < n_bw; ++j){
        auto it_elems = it_elems_begin + j;
        std::vector<smoother_output::pair> &transition_pairs =
          it_elems->transition_pairs;

        it_elems->p = get_cloud_idx(sm_cloud, *(it_bw_begin + j));
        it_elems->log_weight = it_elems->p->log_weight;

        transition_pairs.resize(n_fw);
        auto it_trans_pair = transition_pairs.begin();
        // Loop over second index of pair of particles
        const int *it_fw = it_fw_begin + j * n_fw;
        const double *it_w = it_w_begin + j * n_fw;
        for(unsigned int k = 0; k < n_fw; ++k, ++it_fw, ++it_w, ++it_trans_pair){
          it_trans_pair->p = get_cloud_idx(fw_cloud, *it_fw);
          it_trans_pair->log_weight = log(*it_w);
        }
      }
    }
  }

  return ans;
}
