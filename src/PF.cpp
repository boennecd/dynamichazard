#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"

#define BOOT_FILTER "bootstrap_filter"
#define PF_APPROX_CLOUD_MEAN "PF_normal_approx_w_cloud_mean"
#define AUX_APPROX_CLOUD_MEAN "AUX_normal_approx_w_cloud_mean"
#define PF_APPROX_PARTICLE "PF_normal_approx_w_particles"
#define AUX_APPROX_PARTICLE "AUX_normal_approx_w_particles"

#define PF_COMMON_ARGS                                                \
  const int n_fixed_terms_in_state_vec,                               \
  arma::mat &X,                                                       \
  arma::mat &fixed_terms,                                             \
  const arma::vec &tstart,                                            \
  const arma::vec &tstop,                                             \
  const arma::colvec &a_0,                                            \
  const arma::mat &R,                                                 \
  const arma::mat &L,                                                 \
  const arma::vec &m,                                                 \
  arma::mat &Q_0,                                                     \
  arma::mat &Q,                                                       \
  const arma::mat Q_tilde,                                            \
  const Rcpp::List &risk_obj,                                         \
  const arma::mat &F,                                                 \
  const int n_max,                                                    \
  const int n_threads,                                                \
  const arma::vec &fixed_parems,                                      \
  const int N_fw_n_bw,                                                \
  const int N_smooth,                                                 \
  Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold, \
  const int debug,                                                    \
  const int N_first

/* --------------------------------------- */

template<
  template <
    template <bool> class,
    template <bool> class>
  class smoother>
class PF_smooth_smoother_n_dens {
  using bootstrap_filter_sm =
    smoother<
      None_AUX_resampler,
      importance_dens_no_y_dependence>;

  using PF_w_normal_approx_sm =
    smoother<
      None_AUX_resampler,
      importance_dens_normal_approx_w_cloud_mean>;

  using AUX_w_normal_approx_sm =
    smoother<
      AUX_resampler_normal_approx_w_cloud_mean,
      importance_dens_normal_approx_w_cloud_mean>;

  using PF_w_particles_sm =
    smoother<
      None_AUX_resampler,
      importance_dens_normal_approx_w_particles>;

  using AUX_w_particles_sm =
    smoother<
      AUX_resampler_normal_approx_w_particles,
      importance_dens_normal_approx_w_particles>;

public:
  static Rcpp::List compute(
      const PF_data &data, const std::string method,
      pf_base_dens &dens_calc){
    /* Get the smoothed particles at time 1, 2, ..., d */
    smoother_output result;
    if (method == BOOT_FILTER) {
      result = bootstrap_filter_sm::compute(data, dens_calc);

    } else if (method == PF_APPROX_CLOUD_MEAN){
      result = PF_w_normal_approx_sm::compute(data, dens_calc);

    } else if (method == AUX_APPROX_CLOUD_MEAN){
      result = AUX_w_normal_approx_sm::compute(data, dens_calc);

    } else if (method == PF_APPROX_PARTICLE){
      result = PF_w_particles_sm::compute(data, dens_calc);

    }  else if (method == AUX_APPROX_PARTICLE){
      result = AUX_w_particles_sm::compute(data, dens_calc);

    } else {
      std::stringstream stream;
      stream << "method '" << method << "' is not implemented";
      Rcpp::stop(stream.str());
    }

    /* Create output list */
    return(get_rcpp_list_from_cloud(result, &data));
  }
};

class PF_smooth_dens {
  using Fearnhead_O_N  =
    PF_smooth_smoother_n_dens<PF_smoother_Fearnhead_O_N>;

  using Brier_O_N_square  =
    PF_smooth_smoother_n_dens<PF_smoother_Brier_O_N_square>;

public:
  static Rcpp::List compute(
      const PF_data &data, const std::string smoother,
      const std::string method, pf_base_dens &dens_calc){
    Rcpp::List ans;

    if(smoother == "Fearnhead_O_N"){
      ans = Fearnhead_O_N::compute(data, method, dens_calc);

    } else if (smoother == "Brier_O_N_square"){
      ans = Brier_O_N_square::compute(data, method, dens_calc);

    } else {
      std::stringstream stream;
      stream << "smoother '" << smoother << "' is not implemented";
      Rcpp::stop(stream.str());

    }

    return ans;
  }
};


// [[Rcpp::export]]
Rcpp::List PF_smooth(
    const int n_fixed_terms_in_state_vec, arma::mat &X, arma::mat &fixed_terms,
    const arma::vec &tstart, const arma::vec &tstop, const arma::colvec &a_0,
    const arma::mat &R, const arma::mat &L, const arma::vec &m, arma::mat &Q_0,
    arma::mat &Q, const arma::mat Q_tilde, const Rcpp::List &risk_obj,
    const arma::mat &F, const int n_max, const int n_threads,
    const arma::vec &fixed_parems, const int N_fw_n_bw, const int N_smooth,
    Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
    const int debug, const int N_first,

    /* non-common arguments */
    const std::string method, const std::string smoother, const std::string model){
  const arma::ivec is_event_in_bin = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  random_walk<PF_data> data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin, a_0, R, L, m, Q_0, Q,
      risk_obj, F, n_max, n_threads, fixed_parems, Q_tilde, N_fw_n_bw,
      N_smooth, forward_backward_ESS_threshold, debug, N_first);

  Rcpp::List ans;

  std::unique_ptr<pf_base_dens> dens_calc;
  if(model == "logit"){
    dens_calc.reset(new logistic_dens(data));

  } else if (model == "exponential"){
    dens_calc.reset(new exponential_dens(data));

  } else {
    std::stringstream stream;
    stream << "model '" << model << "' is not implemented";
    Rcpp::stop(stream.str());

  }

  ans = PF_smooth_dens::compute(
    data, smoother, method, *dens_calc.get());

  return(ans);
}

/* --------------------------------------- */


template<bool is_forward>
class PF_single_direction {
  using bootstrap_filter =
    AUX_PF<
      None_AUX_resampler,
      importance_dens_no_y_dependence,
      is_forward>;

  using PF_w_normal_approx =
    AUX_PF<
      None_AUX_resampler,
      importance_dens_normal_approx_w_cloud_mean,
      is_forward>;

  using AUX_w_normal_approx =
    AUX_PF<
      AUX_resampler_normal_approx_w_cloud_mean,
      importance_dens_normal_approx_w_cloud_mean,
      is_forward>;

  using PF_w_particles =
    AUX_PF<
      None_AUX_resampler,
      importance_dens_normal_approx_w_particles,
      is_forward>;

  using AUX_w_particles =
    AUX_PF<
      AUX_resampler_normal_approx_w_particles,
      importance_dens_normal_approx_w_particles,
      is_forward>;

public:
  static Rcpp::List compute(
      const PF_data &data, const std::string method,
      pf_base_dens &dens_calc){
    /* Get the smoothed particles at time 1, 2, ..., d */
    std::vector<cloud> result;
    if (method == BOOT_FILTER) {
      result = bootstrap_filter::compute(data, dens_calc);

    } else if (method == PF_APPROX_CLOUD_MEAN){
      result = PF_w_normal_approx::compute(data, dens_calc);

    } else if (method == AUX_APPROX_CLOUD_MEAN){
      result = AUX_w_normal_approx::compute(data, dens_calc);

    } else if (method == PF_APPROX_PARTICLE){
      result = PF_w_particles::compute(data, dens_calc);

    }  else if (method == AUX_APPROX_PARTICLE){
      result = AUX_w_particles::compute(data, dens_calc);

    } else {
      std::stringstream stream;
      stream << "method '" << method << "' is not implemented";
      Rcpp::stop(stream.str());
    }

    /* Create output list */
    return(get_rcpp_list_from_cloud(
        result, !is_forward, data.state_dim, &data));
  }
};

Rcpp::List PF_single_direction_compute(
    const PF_data &data, const bool is_forward, const std::string method,
    pf_base_dens &dens_calc){
  if(is_forward)
    return PF_single_direction<true>::compute(data, method, dens_calc);

  return PF_single_direction<false>::compute(data, method, dens_calc);
}

// [[Rcpp::export]]
Rcpp::List particle_filter(
    const int n_fixed_terms_in_state_vec, arma::mat &X, arma::mat &fixed_terms,
    const arma::vec &tstart, const arma::vec &tstop, const arma::colvec &a_0,
    const arma::mat &R, const arma::mat &L, const arma::vec &m, arma::mat &Q_0,
    arma::mat &Q, const arma::mat Q_tilde, const Rcpp::List &risk_obj,
    const arma::mat &F, const int n_max, const int n_threads,
    const arma::vec &fixed_parems, const int N_fw_n_bw, const int N_smooth,
    Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
    const int debug, const int N_first,

    /* non-common arguments */
    const bool is_forward, const std::string method, const std::string model){
  const arma::ivec is_event_in_bin =
    Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  random_walk<PF_data> data(
      n_fixed_terms_in_state_vec,
      X,
      fixed_terms,
      tstart,
      tstop,
      is_event_in_bin,
      a_0,
      R,
      L,
      m,
      Q_0,
      Q,
      risk_obj,
      F,
      n_max,
      n_threads,
      fixed_parems,
      Q_tilde,
      N_fw_n_bw,
      N_smooth,
      forward_backward_ESS_threshold,
      debug,
      N_first);

  Rcpp::List ans;

  std::unique_ptr<pf_base_dens> dens_calc;
  if(model == "logit"){
    dens_calc.reset(new logistic_dens(data));

  } else if (model == "exponential"){
    dens_calc.reset(new exponential_dens(data));

  } else {
    std::stringstream stream;
    stream << "model '" << model << "' is not implemented";
    Rcpp::stop(stream.str());

  }

  ans =
    PF_single_direction_compute(data, is_forward, method, *dens_calc.get());

  return(ans);
}

// [[Rcpp::export]]
Rcpp::List compute_summary_stats_first_o_RW(
    const Rcpp::List &rcpp_list, unsigned int n_threads,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0){
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  PF_summary_stats_RW stats;
  {
    auto sm_output = get_clouds_from_rcpp_list(rcpp_list);

    stats = compute_summary_stats_first_o_RW(sm_output, a_0, Q, Q_0);
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
