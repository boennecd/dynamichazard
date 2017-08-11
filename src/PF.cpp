#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"


using PF_easy_fw = AUX_PF<None_AUX_resampler, importance_dens_no_y_dependence, binary, true>;
using PF_easy_bw = AUX_PF<None_AUX_resampler, importance_dens_no_y_dependence, binary, false>;
using PF_easy_sm = PF_smoother<None_AUX_resampler, importance_dens_no_y_dependence, binary>;

using AUX_scale_weights_fw = AUX_PF<AUX_resampler, importance_dens_no_y_dependence, binary, true>;
using AUX_scale_weights_bw = AUX_PF<AUX_resampler, importance_dens_no_y_dependence, binary, false>;
using AUX_scale_weights_sm = PF_smoother<AUX_resampler, importance_dens_no_y_dependence, binary>;

using AUX_scale_weights_n_normal_approx_fw =
  AUX_PF<AUX_resampler, importance_dens_normal_approx, binary, true>;
using AUX_scale_weights_n_normal_approx_bw =
  AUX_PF<AUX_resampler, importance_dens_normal_approx, binary, false>;
using AUX_scale_weights_n_normal_approx_sm =
  PF_smoother<AUX_resampler, importance_dens_normal_approx, binary>;

/* Function to turn a clouds of particles into an Rcpp::List */
Rcpp::List get_rcpp_out_list(
    const PF_data &data, std::vector<cloud> clouds, const bool reverse){
  if(data.debug > 0)
    data.log(1) << "Creating output list";

  auto n_clouds = clouds.size();
  Rcpp::List ans(n_clouds);
  std::vector<cloud>::iterator it =
    reverse ? clouds.end() - 1 : clouds.begin();
  for(std::vector<cloud>::size_type j = 0;
      j < n_clouds;
      ++j, it += 1 - 2 * reverse){
    auto n_elem = it->size();
    arma::uvec parent_idx(n_elem);
    arma::uvec child_idx(n_elem);
    arma::vec weights(n_elem);
    arma::mat states(data.a_0.n_elem, n_elem);

    auto idx_pr = parent_idx.begin();
    auto idx_ch = child_idx.begin();
    auto w = weights.begin();
    auto pr = it->begin();
    for(arma::uword i = 0; i < n_elem; ++i, ++w, ++idx_pr, ++idx_ch, ++pr){
      *idx_pr = (pr->parent) ?
        pr->parent->cloud_idx + 1 /* want non-zero based */ : 0;

      *idx_ch = (pr->child) ?
        pr->child->cloud_idx + 1 /* want non-zero based */ : 0;

      *w = exp(pr->log_weight);
      states.col(i) = pr->state;
    }

    ans[j] =
      Rcpp::List::create(
        Rcpp::Named("parent_idx") = Rcpp::wrap(parent_idx),
        Rcpp::Named("child_idx") = Rcpp::wrap(child_idx),
        Rcpp::Named("weights") = Rcpp::wrap(weights),
        Rcpp::Named("states") = Rcpp::wrap(states)
      );
  }

  return ans;
}

// TODO: want to export?
// [[Rcpp::export]]
Rcpp::List FW_filter(
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
    const int N_first){
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
  std::vector<cloud> clouds = AUX_scale_weights_fw::compute(data);

  /* Create output list */
  return(get_rcpp_out_list(data, clouds, false));
}

// TODO: want to export?
// [[Rcpp::export]]
Rcpp::List BW_filter(
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
    const int N_first){
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
  std::vector<cloud> clouds = AUX_scale_weights_bw::compute(data);

  /* Create output list */
  return(get_rcpp_out_list(data, clouds, true));
}

// TODO: want to export?
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
  if(method == "AUX"){
    result = AUX_scale_weights_sm::compute(data);
  } else if (method == "PF") {
    result = PF_easy_sm::compute(data);
  } else if (method == "normal_approx"){
    result = AUX_scale_weights_n_normal_approx_sm::compute(data);
  }else
    Rcpp::stop("'method' not implemented");


  /* Create output list */
  return Rcpp::List::create(
    Rcpp::Named("forward_clouds") =
      get_rcpp_out_list(data, result.forward_clouds, false),
    Rcpp::Named("backward_clouds") =
      get_rcpp_out_list(data, result.backward_clouds, true),
    Rcpp::Named("smoothed_clouds") =
      get_rcpp_out_list(data, result.smoothed_clouds, false));
}
