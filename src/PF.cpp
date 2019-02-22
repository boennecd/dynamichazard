#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"
#include "PF/est_params.h"

#define BOOT_FILTER "bootstrap_filter"
#define PF_APPROX_CLOUD_MEAN "PF_normal_approx_w_cloud_mean"
#define AUX_APPROX_CLOUD_MEAN "AUX_normal_approx_w_cloud_mean"
#define PF_APPROX_PARTICLE "PF_normal_approx_w_particles"
#define AUX_APPROX_PARTICLE "AUX_normal_approx_w_particles"

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
      pf_dens &dens_calc){
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

std::string get_family(const std::string &model){
  if(model == "logit"){
    return BINOMIAL;

  } else if (model == "exponential"){
    return POISSON;

  } else if (model == "cloglog") {
    return CLOGLOG;

  }

  std::stringstream stream;
  stream << "model '" << model << "' is not implemented";
  Rcpp::stop(stream.str());

  return "";
}

class PF_smooth_dens {
  using Fearnhead_O_N  =
    PF_smooth_smoother_n_dens<PF_smoother_Fearnhead_O_N>;

  using Brier_O_N_square  =
    PF_smooth_smoother_n_dens<PF_smoother_Brier_O_N_square>;

public:
  static Rcpp::List compute(
      const PF_data &data, const std::string smoother,
      const std::string method, pf_dens &dens_calc){
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
    const arma::mat &R, arma::mat &Q_0,
    arma::mat &Q, const arma::mat Q_tilde, const Rcpp::List &risk_obj,
    const arma::mat &F, const int n_max, const int n_threads,
    const arma::vec &fixed_params, const int N_fw_n_bw, const int N_smooth,
    const int N_smooth_final, const double covar_fac, const double ftol_rel,
    Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
    const int debug, const int N_first, std::string type, const int nu,

    /* non-common arguments */
    const std::string method, const std::string smoother, const std::string model){
  const arma::ivec is_event_in_bin = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  std::unique_ptr<PF_data> data(new PF_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin, a_0, R, R.t(), Q_0, Q,
      risk_obj, F, n_max, n_threads, fixed_params, Q_tilde,
      N_fw_n_bw, N_smooth, N_smooth_final, forward_backward_ESS_threshold,
      debug, N_first, nu, covar_fac, ftol_rel));

  Rcpp::List ans;
  std::string family_use = get_family(model);
  pf_dens dens_calc(*data, family_use);
  return PF_smooth_dens::compute(
    *data.get(), smoother, method, dens_calc);
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
      pf_dens &dens_calc){
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
    pf_dens &dens_calc){
  if(is_forward)
    return PF_single_direction<true>::compute(data, method, dens_calc);

  return PF_single_direction<false>::compute(data, method, dens_calc);
}

// [[Rcpp::export]]
Rcpp::List particle_filter(
    const int n_fixed_terms_in_state_vec, arma::mat &X, arma::mat &fixed_terms,
    const arma::vec &tstart, const arma::vec &tstop, const arma::colvec &a_0,
    const arma::mat &R, arma::mat &Q_0,
    arma::mat &Q, const arma::mat Q_tilde, const Rcpp::List &risk_obj,
    const arma::mat &F, const int n_threads,
    const arma::vec &fixed_params, const int N_fw_n_bw,
    Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
    const int debug, const int N_first, const int nu, std::string type,
    const bool is_forward, const std::string method, const std::string model,
    const double covar_fac, const double ftol_rel){
  const arma::ivec is_event_in_bin =
    Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  const unsigned int N_smooth = 1, n_max = 1, N_smooth_final = 1;
  std::unique_ptr<PF_data> data(new PF_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin, a_0, R, R.t(), Q_0, Q,
      risk_obj, F, n_max, n_threads, fixed_params, Q_tilde,
      N_fw_n_bw, N_smooth, N_smooth_final, forward_backward_ESS_threshold,
      debug, N_first, nu, covar_fac, ftol_rel));

  std::string family_use = get_family(model);
  pf_dens dens_calc(*data, family_use);
  return PF_single_direction_compute(
    *data.get(), is_forward, method, dens_calc);
}

// [[Rcpp::export]]
Rcpp::List compute_PF_summary_stats(
    const Rcpp::List &rcpp_list, unsigned int n_threads,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
    const arma::mat &R, const bool debug, const arma::mat F,
    const bool do_use_F = false, const bool do_compute_E_x = true){
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  PF_summary_stats stats;
  {
    auto sm_output = get_clouds_from_rcpp_list(rcpp_list);

    stats = compute_PF_summary_stats(
      sm_output, a_0, Q, Q_0, F, do_use_F, do_compute_E_x);
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

// [[Rcpp::export]]
Rcpp::List PF_est_params_dens(
    const Rcpp::List &rcpp_list, unsigned int n_threads,
    const arma::vec &a_0, const arma::mat &Q, const arma::mat &Q_0,
    const arma::mat &R, const bool debug, const bool do_est_a_0 = false,
    const bool only_QR = false){
  const unsigned long int max_bytes = 5000000;
  PF_parameters new_params;

  {
    auto sm_output = get_clouds_from_rcpp_list(rcpp_list);

    new_params = est_params_dens(
      sm_output, a_0, Q, Q_0, R, n_threads, do_est_a_0, debug, max_bytes,
      only_QR);
  }

  return Rcpp::List::create(
    Rcpp::Named("a_0")     = new_params.a_0,
    Rcpp::Named("R_top_F") = new_params.R_top_F,
    Rcpp::Named("Q")       = new_params.Q,
    Rcpp::Named("QR_R")    = new_params.R,
    Rcpp::Named("QR_F")    = new_params.F,
    Rcpp::Named("QR_dev")  = new_params.dev);
}



// [[Rcpp::export]]
Rcpp::List PF_get_score_n_hess_cpp(
  const Rcpp::List fw_cloud, const arma::mat &Q,
  const arma::mat &F, Rcpp::List risk_obj,
  const arma::mat &ran_vars, const arma::mat &fixed_terms,
  const arma::vec &tstart, const arma::vec &tstop,
  const arma::vec &fixed_params, const std::string family,
  const int max_threads, const bool debug, const bool only_score = false)
{
  std::vector<cloud> clouds_cpp =
    get_cloud_from_rcpp_list<false, false>(fw_cloud);
  const Rcpp::List risk_sets_R = Rcpp::as<Rcpp::List>(risk_obj["risk_sets"]);
  std::vector<arma::uvec> risk_sets(risk_sets_R.size());
  for(unsigned int i = 0; i < risk_sets_R.size(); ++i)
    risk_sets[i] = Rcpp::as<arma::uvec>(risk_sets_R[i]) - 1L;

  arma::ivec is_event_in = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);
  arma::vec event_times = Rcpp::as<arma::vec>(risk_obj["event_times"]);

  std::vector<score_n_hess> out_cpp = PF_get_score_n_hess
    (clouds_cpp, Q, F, risk_sets, is_event_in, event_times, ran_vars,
     fixed_terms, tstart, tstop, fixed_params, family, max_threads, debug,
     only_score);

  arma::uword
    p = out_cpp[0].get_a_state().n_elem, q = out_cpp[0].get_a_obs().n_elem,
    n = out_cpp.size();
  int i_p = p, i_q = q;
  arma::vec S_state(p, arma::fill::zeros), S_obs(q, arma::fill::zeros);
  arma::mat neg_obs_info_state_dd(p, p, arma::fill::zeros),
            neg_obs_info_state_d (p, p, arma::fill::zeros),
            neg_obs_info_obs_dd  (q, q, arma::fill::zeros),
            neg_obs_info_obs_d   (q, q, arma::fill::zeros);

  /* TODO: issue w/ catastrophic cancellation? */
  auto pi = clouds_cpp.back().begin();
  for(auto o = out_cpp.begin(); o != out_cpp.end(); ++o, ++pi){
    double w = exp(pi->log_weight);
    S_state += w * o->get_a_state();
    S_obs   += w * o->get_a_obs();

    if(!only_score){
      R_BLAS_LAPACK::sym_mat_rank_one_update
        (&i_p, &w, o->get_a_state().memptr(), neg_obs_info_state_d.memptr());
      R_BLAS_LAPACK::sym_mat_rank_one_update
        (&i_q, &w, o->get_a_obs()  .memptr(), neg_obs_info_obs_d  .memptr());

      neg_obs_info_state_dd += w * o->get_B_state();
      neg_obs_info_obs_dd   += w * o->get_B_obs();
    }
  }

  if(!only_score){
    neg_obs_info_state_d    = arma::symmatu(neg_obs_info_state_d);
    neg_obs_info_state_dd   = arma::symmatu(neg_obs_info_state_dd);
    neg_obs_info_obs_d      = arma::symmatu(neg_obs_info_obs_d);
    neg_obs_info_obs_dd     = arma::symmatu(neg_obs_info_obs_dd);

  } else {
    neg_obs_info_state_d.fill(NA_REAL);
    neg_obs_info_state_dd.fill(NA_REAL);
    neg_obs_info_obs_d.fill(NA_REAL);
    neg_obs_info_obs_dd.fill(NA_REAL);

  }

  return Rcpp::List::create(
    Rcpp::Named("S_state")              = std::move(S_state),
    Rcpp::Named("S_obs")                = std::move(S_obs),
    Rcpp::Named("E_second_deriv_state") = std::move(neg_obs_info_state_dd),
    Rcpp::Named("E_second_deriv_obs")   = std::move(neg_obs_info_obs_dd),
    Rcpp::Named("E_score_outer_state")  = std::move(neg_obs_info_state_d),
    Rcpp::Named("E_score_outer_obs")    = std::move(neg_obs_info_obs_d));
}
