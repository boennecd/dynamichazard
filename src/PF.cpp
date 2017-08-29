#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"

using bootstrap_filter_sm = PF_smoother<
  None_AUX_resampler,
  importance_dens_no_y_dependence,
  binary>;

using PF_w_normal_approx_sm =
  PF_smoother<
    None_AUX_resampler,
    importance_dens_normal_approx_w_cloud_mean,
    binary>;

using AUX_w_normal_approx_sm =
  PF_smoother<
    AUX_resampler_normal_approx_w_cloud_mean,
    importance_dens_normal_approx_w_cloud_mean,
    binary>;

using PF_w_particles_sm =
  PF_smoother<
    None_AUX_resampler,
    importance_dens_normal_approx_w_particles,
    binary>;

using AUX_w_particles_sm =
  PF_smoother<
    AUX_resampler_normal_approx_w_particles,
    importance_dens_normal_approx_w_particles,
    binary>;

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
      *idx_pr = (pr->parent) ?
        pr->parent->cloud_idx + 1 /* want non-zero based */ : 0;

      *idx_ch = (pr->child) ?
        pr->child->cloud_idx + 1 /* want non-zero based */ : 0;

      *w = exp(pr->log_weight);
      *un_w = pr->log_unnormalized_weight;

      states.col(i) = pr->state;
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
  unsigned int state_dim = sm_output.forward_clouds[0][0].state.n_elem;

  Rcpp::List ans = Rcpp::List::create(
    Rcpp::Named("forward_clouds") =
      get_rcpp_list_from_cloud(
        sm_output.forward_clouds, false, state_dim, data),

    Rcpp::Named("backward_clouds") =
      get_rcpp_list_from_cloud(
        sm_output.backward_clouds, true, state_dim, data),

    Rcpp::Named("smoothed_clouds") =
      get_rcpp_list_from_cloud(
        sm_output.smoothed_clouds, false, state_dim, data));

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
          parent = &(*fw)[i /* starts at time 0 */][*it_par - 1];

        }else
          parent = &(*it_ans_prev)[*it_par - 1];

        const particle *child  =
          /* starts at time 1 and has on more index */
          is_smooth ? &(*bw)[n_periods - (i + 1)][*it_child - 1] : nullptr;

        particle &p = it_ans->New_particle(states.col(j), parent, child);
        p.log_weight = log(*it_w);
        p.log_unnormalized_weight = *it_un_w;
      }

      if(reverse && i == 0) break;
    }

    return ans;
}

smoother_output get_clouds_from_rcpp_list(const Rcpp::List &rcpp_list){
  smoother_output ans;

  std::vector<cloud> &fw_clouds = ans.forward_clouds;
  std::vector<cloud> &bw_clouds = ans.backward_clouds;
  std::vector<cloud> &smooth_clouds = ans.smoothed_clouds;

  fw_clouds = get_clouds_from_rcpp_list_util<false, false>(rcpp_list["forward_clouds"]);
  bw_clouds = get_clouds_from_rcpp_list_util<false, true>(rcpp_list["backward_clouds"]);
  smooth_clouds = get_clouds_from_rcpp_list_util<true, false>(
    rcpp_list["smoothed_clouds"], &fw_clouds, &bw_clouds);

  return ans;
}

/* --------------------------------------- */

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
    const std::string method){
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

/* --------------------------------------- */

struct PF_summary_stats {
  /* E[X_t] */
  std::vector<arma::vec> E_xs;
  /* E[(X_t - FX_{t-1})(X_t - FX_{t-1})^T] */
  std::vector<arma::mat> E_x_less_x_less_one_outers;
};

static PF_summary_stats compute_summary_stats(const std::vector<cloud> &smoothed_clouds){
  PF_summary_stats ans;
  std::vector<arma::vec> &E_xs = ans.E_xs;
  std::vector<arma::mat> &E_x_less_x_less_one_outers = ans.E_x_less_x_less_one_outers;

  unsigned int n_periods = smoothed_clouds.size();
  unsigned int n_elem = smoothed_clouds[0][0].state.n_elem;

  for(unsigned int i = 0; i < n_periods; ++i){
    const bool is_last = i == n_periods - 1;
    auto it_cl = (smoothed_clouds.begin() + i);

    unsigned int n_part = it_cl->size();
    arma::vec E_x(n_elem, arma::fill::zeros);
    arma::mat E_x_less_x_less_one_outer(n_elem, n_elem, arma::fill::zeros);

    auto cl_begin = it_cl->begin();
    for(unsigned int j = 0; j < n_part; ++j){
      const particle &p = *(cl_begin + j);
      double weight = exp(p.log_weight);

      E_x += weight * p.state;
      // TODO: use BLAS outer product rank one update
      if(is_last)
        continue;
      arma::vec inter = sqrt(weight) * (p.state - p.parent->state);
      E_x_less_x_less_one_outer += inter * inter.t();
    }

    if(is_last){
      --it_cl;
      for(auto p = it_cl->begin(); p != it_cl->end(); ++p){
        // TODO: use BLAS outer product rank one update
        arma::vec inter =
          sqrt(exp(p->log_weight)) * (p->child->state - p->state);
        E_x_less_x_less_one_outer += inter * inter.t();
      }
    }

    E_xs.push_back(std::move(E_x));
    E_x_less_x_less_one_outers.push_back(std::move(E_x_less_x_less_one_outer));
  }

  return ans;
}

static PF_summary_stats compute_summary_stats(const smoother_output sm_output){
  return(compute_summary_stats(sm_output.smoothed_clouds));
}

// [[Rcpp::export]]
Rcpp::List compute_summary_stats(const Rcpp::List &rcpp_list){
  PF_summary_stats stats;
  {
    auto sm_output = get_clouds_from_rcpp_list(rcpp_list);
    stats = compute_summary_stats(sm_output);
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


