#include "PF/PFs.h"
#include "PF/importance_samplers.h"
#include "PF/resamplers.h"
#include "PF/densities.h"

using simple_smoother = PF_smoother<None_AUX_resampler, importance_dens_no_y_dependence, binary>;

// TODO: want to export?
// [[Rcpp::export]]
Rcpp::List PF_smooth(
    const int n_fixed_terms_in_state_vec_,
    arma::mat &X_,
    arma::mat &fixed_terms_,
    const arma::vec &tstart_,
    const arma::vec &tstop_, const arma::ivec &is_event_in_bin_,
    const arma::colvec &a_0,
    arma::mat &Q_0_,
    arma::mat &Q_,
    const Rcpp::List &risk_obj,
    const arma::mat &F__,
    const int n_max,
    const int order_,
    const int n_threads_,
    const int N_fw_n_bw,
    const int N_smooth,
    const double forward_backward_ESS_threshold){

  PF_data data(
      n_fixed_terms_in_state_vec_,
      X_,
      fixed_terms_,
      tstart_,
      tstop_, is_event_in_bin_,
      a_0,
      Q_0_,
      Q_,
      risk_obj,
      F__,
      n_max,
      order_,
      n_threads_,

      // new arguments
      N_fw_n_bw,
      N_smooth,
      forward_backward_ESS_threshold);

  std::vector<cloud> smoothed_cloud = simple_smoother::compute(data);
}

// exported for test
// [[Rcpp::export]]
double dmvnrm_log_test(
    const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv){
  return(dmvnrm_log(x, mean, sigma_chol_inv));
}
