#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"
#include "PF/est_params.h"

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
  const Rcpp::List fw_cloud, arma::mat &Q,
  const arma::mat &F, Rcpp::List risk_obj,
  arma::mat &ran_vars, arma::mat &fixed_terms,
  const arma::vec &tstart, const arma::vec &tstop,
  const arma::vec &fixed_params, const std::string family,
  const int max_threads, const bool debug, const arma::vec &a_0,
  const arma::mat &R, arma::mat &Q_0, const arma::mat &Q_tilde,
  const arma::uword N_fw_n_bw, const arma::uword N_first, const double nu,
  const double covar_fac, const double ftol_rel, const std::string method,
  Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
  const bool use_O_n_sq = false, const bool only_score = false)
{
  std::vector<cloud> clouds_cpp;
  const Rcpp::List risk_sets_R = Rcpp::as<Rcpp::List>(risk_obj["risk_sets"]);
  std::vector<arma::uvec> risk_sets(risk_sets_R.size());
  for(unsigned int i = 0; i < risk_sets_R.size(); ++i)
    risk_sets[i] = Rcpp::as<arma::uvec>(risk_sets_R[i]) - 1L;

  arma::ivec is_event_in = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);
  arma::vec event_times = Rcpp::as<arma::vec>(risk_obj["event_times"]);

  /* get intermediary values */
  std::vector<std::unique_ptr<score_n_hess_base> > out_cpp;
  if(use_O_n_sq){
    out_cpp = PF_get_score_n_hess_O_N_sq
    (Q, F, risk_sets, risk_obj, is_event_in, event_times, ran_vars,
     fixed_terms, tstart, tstop, fixed_params, family, max_threads, debug,
     only_score, a_0, R, Q_0, Q_tilde, N_fw_n_bw, N_first, nu, covar_fac,
     ftol_rel, forward_backward_ESS_threshold, method);

  } else {
    clouds_cpp = get_cloud_from_rcpp_list<false, false>(fw_cloud);
    out_cpp = PF_get_score_n_hess
      (clouds_cpp, Q, F, risk_sets, is_event_in, event_times, ran_vars,
       fixed_terms, tstart, tstop, fixed_params, family, max_threads, debug,
       only_score);
  }

  /* compute score vector, weigthed sum of outer products of scores, and
   * weighted sum of hessian terms */
  const int dscore = out_cpp[0]->get_score().n_elem;
  arma::vec score      (dscore, arma::fill::zeros);
  arma::mat score_outer(dscore, dscore, arma::fill::zeros),
            hess_terms (dscore, dscore, arma::fill::zeros);

  for(auto &o : out_cpp){
    const double w = o->get_weight();
    score += w * o->get_score();

    if(!only_score){
      R_BLAS_LAPACK::sym_mat_rank_one_update
        (&dscore, &w, o->get_score().memptr(), score_outer.memptr());
      hess_terms += w * o->get_hess_terms();

    }
  }

  if(!only_score){
    score_outer    = arma::symmatu(score_outer);
    hess_terms     = arma::symmatu(hess_terms);

  } else {
    score_outer.fill(NA_REAL);
    hess_terms.fill(NA_REAL);

  }

  return Rcpp::List::create(
    Rcpp::Named("score")       = std::move(score),
    Rcpp::Named("score_outer") = std::move(score_outer),
    Rcpp::Named("hess_terms")  = std::move(hess_terms));
}
