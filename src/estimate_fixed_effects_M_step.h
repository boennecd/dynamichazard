#ifndef DDHAZARD_EST_M_STEP_H
#define DDHAZARD_EST_M_STEP_H
#include "bigglm_wrapper.h"
#include "problem_data.h"
#include "arma_n_rcpp.h"

// Method to estimate fixed effects like in biglm::bigglm
template <class updater>
void estimate_fixed_effects_M_step(ddhazard_data * const p_data, arma::uword chunk_size){
  using uword = arma::uword;

  // Data storage class
  struct obs_info {
    uword data_row_idx;
    bool is_event_in_bin;
    double offset;
    double weight;
  };

  // Set up loop variables
  std::vector<obs_info> obs_data;
  double bin_stop = p_data->min_start;
  int t = 1;
  for(auto it = p_data->risk_sets.begin();
      it != p_data->risk_sets.end(); ++it, ++t){
    double bin_start = bin_stop;
    double delta_t = p_data->I_len[t - 1]; // I_len is std::vector and thus uses zero index
    bin_stop += delta_t;

    // Find the risk set
    arma::uvec r_set = Rcpp::as<arma::uvec>(*it) - 1;

    // Compute offsets from dynamic effects if needed
    arma::vec offsets;
    if(p_data->any_dynamic){
      offsets =
        // .col(t) and not .col(t - 1) this is a_(t - 1 | d)
        p_data->X.cols(r_set).t() * p_data->a_t_t_s.col(t).head(p_data->n_params_state_vec);
    } else {
      offsets = arma::vec(r_set.n_elem, arma::fill::zeros);
    }

    // Add elements to vector
    auto off_i = offsets.begin();
    for(auto i = r_set.begin(); i != r_set.end(); ++i, ++off_i){
      bool is_event_in_bin  = p_data->is_event_in_bin(*i) == (t - 1);
      double offset =
        *off_i + updater::family::time_offset(
            std::min(p_data->tstop(*i), bin_stop) -
            std::max(p_data->tstart(*i), bin_start));
      double w = p_data->weights(*i);
      obs_info obj = {*i, is_event_in_bin, offset, w};
      obs_data.push_back(std::move(obj));
    }
  }

  /* Tried to use sampling because there as I suspected it might cause a
   * I guess there is not an error (or just a small one) of the series of rank
   * one updates in the QR decomposition. I have kept the code and comments for
   * now...
   *
   * I should likely have looked at source like 6.5 of: Golub, Gene H. Matrix Computations (Johns Hopkins Studies in the Mathematical Sciences) (Kindle Location 8376). Johns Hopkins University Press. Kindle Edition.
   *
   * Sample obs_info. This is done to avoid issues if one of the fixed
   * covariates depends on time too.
   * I use the sampler from http://gallery.rcpp.org/articles/using-the-Rcpp-based-sample-implementation/
   * The motivation is to get reproducible results. See the post here for
   * using the random number generator with Rcpp: http://gallery.rcpp.org/articles/random-number-generation/

   arma::uvec ord = arma::linspace<arma::uvec>(0, obs_data.size() - 1L, obs_data.size());
   ord = Rcpp::RcppArmadillo::sample(ord, ord.n_elem, false);

   */

  // Start the estimation
  int it_outer = 0;
  arma::vec old_beta;
  uword n_elements = obs_data.size();
  do{
    qr_obj qr(p_data->fixed_parems.n_elem);

    for(uword i = 0; i < n_elements; ){
      // Find number of elements in this block and define initialize obects
      uword n_take = std::min<uword>(chunk_size, n_elements - i);
      arma::vec offsets(n_take);
      arma::mat fixed_terms(p_data->fixed_parems.n_elem, n_take);
      arma::vec y(n_take);
      arma::vec w(n_take);

      auto o_j = offsets.begin();
      auto y_j = y.begin();
      auto w_j = w.begin();
      for(uword j = 0; j < n_take; ++j, ++i, ++o_j, ++y_j, ++w_j){
        obs_info *dat = &obs_data[i];
        *o_j = dat->offset;
        *y_j = dat->is_event_in_bin;
        *w_j = dat->weight;
        fixed_terms.col(j) = p_data->fixed_terms.col(dat->data_row_idx);
      }
      arma::vec eta = fixed_terms.t() * p_data->fixed_parems;

      // Update QR decomposition
      updater::update(qr, fixed_terms, eta, offsets, y, w);
    }

    old_beta = p_data->fixed_parems;
    p_data->fixed_parems = bigglm_regcf(qr);

    if(p_data->debug){
      my_debug_logger(*p_data) << "Iteration " << it_outer + 1 << " of estimating fixed effects in M-step";
      my_print(*p_data, old_beta, "Fixed effects before update");
      my_print(*p_data, p_data->fixed_parems, "Fixed effects after update");
    }

  } while(++it_outer < p_data->max_it_fixed_params && // Key that this the first condition we check when we use &&
    arma::norm(p_data->fixed_parems - old_beta, 2) / (arma::norm(old_beta, 2) + 1e-8) > p_data->eps_fixed_parems);

  static bool failed_to_converge_once = false;
  if(it_outer == p_data->max_it_fixed_params && !failed_to_converge_once){
    failed_to_converge_once = true;
    std::stringstream msg;
    msg << "Failed to estimate fixed effects in " << p_data->max_it_fixed_params << " iterations at least once" << std::endl;
    Rcpp::warning(msg.str());
  }
}
#endif
