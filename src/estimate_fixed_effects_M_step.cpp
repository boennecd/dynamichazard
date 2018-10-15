#include "estimate_fixed_effects_M_step.h"
#include "bigglm_wrapper.h"
#include "utils.h"

void estimate_fixed_effects_M_step(
    ddhazard_data * const p_data, arma::uword chunk_size, family_base &fam){
  using uword = arma::uword;

  // Data storage class
  struct obs_info {
    uword data_row_idx;
    bool is_event_in_bin;
    double offset;
    double at_risk_length;
    double weight;
  };

  // Set up loop variables
  std::vector<obs_info> obs_data;
  double bin_stop = p_data->min_start;
  int t = 1;
  const bool uses_at_risk_length = fam.uses_at_risk_length();
  for(auto it = p_data->risk_sets.begin();
      it != p_data->risk_sets.end(); ++it, ++t){
    double bin_start = bin_stop;
    double delta_t = p_data->I_len[t - 1]; // I_len is std::vector and thus uses zero index
    bin_stop += delta_t;

    // Find the risk set
    arma::uvec r_set = get_risk_set(*p_data, t);

    // Compute offsets from dynamic effects if needed
    arma::vec offsets;
    if(p_data->any_dynamic){
      offsets =
        p_data->X.cols(r_set).t() *
        p_data->state_lp->map(p_data->a_t_t_s.col(t)).sv;

    } else {
      offsets = arma::vec(r_set.n_elem, arma::fill::zeros);

    }

    // Add elements to vector
    auto off_i = offsets.begin();
    for(auto i = r_set.begin(); i != r_set.end(); ++i, ++off_i){
      bool is_event_in_bin  = p_data->is_event_in_bin(*i) == (t - 1);
      double offset = *off_i;

      double at_risk_length =
        uses_at_risk_length ?
        get_at_risk_length(
          p_data->tstop(*i), bin_stop, p_data->tstart(*i), bin_start) : 0;
      double w = p_data->weights(*i);
      obs_info obj = {*i, is_event_in_bin, offset, at_risk_length, w};
      obs_data.push_back(std::move(obj));
    }
  }

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
      arma::vec at_risk_length(n_take);
      arma::vec y(n_take);
      arma::vec w(n_take);

      auto o_j = offsets.begin();
      double *at_j = at_risk_length.memptr();
      auto y_j = y.begin();
      auto w_j = w.begin();
      for(uword j = 0; j < n_take; ++j, ++i, ++o_j, ++at_j, ++y_j, ++w_j){
        obs_info *dat = &obs_data[i];
        *o_j  = dat->offset;
        *at_j = dat->at_risk_length;
        *y_j  = dat->is_event_in_bin;
        *w_j  = dat->weight;
        fixed_terms.col(j) = p_data->fixed_terms.col(dat->data_row_idx);
      }
      arma::vec eta = fixed_terms.t() * p_data->fixed_parems;

      // Update QR decomposition
      bigglm_updateQR::update(
        qr, fixed_terms, eta, offsets, at_risk_length, y, w, fam);
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
