#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"

using simple_smoother = PF_smoother<None_AUX_resampler, importance_dens_no_y_dependence, binary>;

/* Function to turn a clouds of particles into an Rcpp::List */
Rcpp::List get_rcpp_out_list(const PF_data &data, std::vector<cloud> clouds){
  if(data.debug > 0)
    data.log(1) << "Creating output list";

  auto n_clouds = clouds.size();
  Rcpp::List ans(n_clouds);
  auto it = clouds.begin();
  for(std::vector<cloud>::size_type j = 0; j < n_clouds; ++j, ++it){
    auto n_elem = it->size();
    arma::uvec parent_idx(n_elem);
    arma::vec weights(n_elem);
    arma::mat states(data.a_0.n_elem, n_elem);

    auto idx = parent_idx.begin();
    auto w = weights.begin();
    auto pr = it->begin();
    for(arma::uword i = 0; i < n_elem; ++i, ++idx, ++w, ++pr){

      if(j > 0) {
        *idx = pr->parent->cloud_idx + 1; /* want non-zero based */
      } else {
        *idx = 0;
      }

      *w = exp(pr->log_weight);
      states.col(i) = pr->state;
    }

    ans[j] =
      Rcpp::List::create(
        Rcpp::Named("parent_idx") = Rcpp::wrap(parent_idx),
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
    const Rcpp::List &risk_obj,
    const arma::mat &F,
    const int n_max,
    const int order,
    const int n_threads,
    const int N_fw_n_bw,
    const int N_smooth,
    Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
    const int debug){
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
      N_fw_n_bw,
      N_smooth,
      forward_backward_ESS_threshold,
      debug);

  /* Get the smoothed particles at time 1, 2, ..., d */
  std::vector<cloud> clouds = AUX_PF<
    None_AUX_resampler, importance_dens_no_y_dependence, binary, true>::compute(data);

  /* Create output list */
  return(get_rcpp_out_list(data, clouds));
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
    const Rcpp::List &risk_obj,
    const arma::mat &F,
    const int n_max,
    const int order,
    const int n_threads,
    const int N_fw_n_bw,
    const int N_smooth,
    Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold,
    const int debug){
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
      N_fw_n_bw,
      N_smooth,
      forward_backward_ESS_threshold,
      debug);

  /* Get the smoothed particles at time 1, 2, ..., d */
  std::vector<cloud> smoothed_clouds = simple_smoother::compute(data);

  /* Create output list */
  return(get_rcpp_out_list(data, smoothed_clouds));
}

// exported for test
// [[Rcpp::export]]
double dmvnrm_log_test(
    const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv){
  return(dmvnrm_log(x, mean, sigma_chol_inv));
}
