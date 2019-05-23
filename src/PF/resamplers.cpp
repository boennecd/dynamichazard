#include "resamplers.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

template <arma::uvec (*sample_func)(const arma::uword, arma::vec&)>
arma::uvec resampler_base<sample_func>::sample(
    const PF_data &data, arma::vec &probs, const double ESS, bool &did_resample){
  if(probs.n_elem != data.N_fw_n_bw){
    if(data.debug > 1){
      data.log(2) << "Subsampling " << probs.n_elem << " to get "
                  << data.N_fw_n_bw <<  " using re-sampling weights. ESS of re-sampling weights are "
                  << ESS;
    }

    did_resample = true;
    return(sample_func(data.N_fw_n_bw, probs));
  }

  if(ESS < data.forward_backward_ESS_threshold){
    if(data.debug > 1){
      data.log(2) << "ESS of re-sampling weights is below threshold (" << ESS << " < "
                  << data.forward_backward_ESS_threshold << "). Re-sampling";
    }

    if(data.debug > 2){
      data.log(3) << "Re-sampling " << data.N_fw_n_bw << " indices "
                  << " from " <<  probs.n_elem << " elements "
                  << " with " <<  arma::max(probs) << " as the higest probability";
    }

    did_resample = true;
    return(sample_func(data.N_fw_n_bw, probs));
  }

  if(data.debug > 1){
    data.log(2) << "ESS of re-sampling weights is greater than threshold (" << ESS << " >= "
                << data.forward_backward_ESS_threshold << "). No re-sampling needed";
  }

  did_resample = false;
  return(arma::linspace<arma::uvec>(0, data.N_fw_n_bw - 1, data.N_fw_n_bw));
}

template class resampler_base<systematic_resampling>;



template<bool is_forward>
nothing None_AUX_resampler<is_forward>::resampler(
    pf_dens &dens_calc, const PF_data &data, cloud &PF_cloud,
    std::shared_ptr<PF_cdist> y_dist,
    arma::uword t, arma::uvec &outcome, bool &did_resample){
  /* Compute effective sample size (ESS) */
  arma::vec weights(PF_cloud.size());
  double ESS = 0;
  auto w = weights.begin();
  for(auto it = PF_cloud.begin(); it != PF_cloud.end(); ++it, ++w){
    /* No need to update weights */
    it->log_resampling_weight = it->log_weight;

    *w = exp(it->log_resampling_weight);
    ESS += *w * *w;
  }
  ESS = 1/ ESS;

  outcome = sample(data, weights, ESS, did_resample);

  return nothing();
}

template class None_AUX_resampler<false>;
template class None_AUX_resampler<true>;



template<bool is_forward>
std::vector<std::unique_ptr<dist_comb>>
AUX_resampler_normal_approx_w_cloud_mean<is_forward>::resampler(
    pf_dens &dens_calc, const PF_data &data, cloud &PF_cloud,
    std::shared_ptr<PF_cdist> y_dist,
    arma::uword t, arma::uvec &outcome, bool &did_resample) {
  get_approx_use_mean_output ans =
    get_approx_use_mean<is_forward>(y_dist, PF_cloud, data, dens_calc, t);

  std::shared_ptr<PF_cdist> prior, prior_p1;
  if(!is_forward){
    prior = dens_calc.get_prior(t);
    prior_p1 = dens_calc.get_prior(t + 1L);
  }

  /* Compute sampling weights */
  double max_weight =  -std::numeric_limits<double>::max();
  unsigned int n_elem = PF_cloud.size();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:max_weight)
#endif
  for(unsigned int i = 0; i < n_elem; ++i){ // loop over cloud elements
    particle &p = PF_cloud[i];
    std::unique_ptr<dist_comb> &dc = ans.dists[i];

    double log_prob_y_given_state = y_dist->log_dens(dc->get_mean()),
      log_prop_transition = is_forward ?
        dens_calc.log_prob_state_given_parent(
          dc->get_mean(), p.get_state()) :
        dens_calc.log_prob_state_given_child(
          dc->get_mean(), p.get_state()),
      log_prop_proposal = dc->log_density(dc->get_mean());

    p.log_resampling_weight =
      p.log_weight + log_prop_transition + log_prob_y_given_state
      - log_prop_proposal;

    if(!is_forward)
      p.log_resampling_weight += (
        prior->log_dens(dc->get_mean()) - prior_p1->log_dens(p.get_state()));

    max_weight = MAX(max_weight, p.log_resampling_weight);
  } // end loop over cloud elements

  auto norm_out =
    normalize_log_resampling_weight
    <true, true>
    (PF_cloud, max_weight);
  outcome = sample(data, norm_out.weights, norm_out.ESS, did_resample);

  return std::move(ans.dists);
}

template class AUX_resampler_normal_approx_w_cloud_mean<false>;
template class AUX_resampler_normal_approx_w_cloud_mean<true>;



template<bool is_forward>
std::vector<std::unique_ptr<dist_comb>>
AUX_resampler_normal_approx_w_particles<is_forward>::resampler(
  pf_dens &dens_calc, const PF_data &data, cloud &PF_cloud,
  std::shared_ptr<PF_cdist> y_dist,
  arma::uword t, arma::uvec &outcome, bool &did_resample) {
  get_approx_use_particle_output ans =
    get_approx_use_particle<is_forward>(y_dist, PF_cloud, data, dens_calc, t);

  std::shared_ptr<PF_cdist> prior, prior_p1;
  if(!is_forward){
    prior = dens_calc.get_prior(t);
    prior_p1 = dens_calc.get_prior(t + 1L);
  }

  /* Compute sampling weights */
  double max_weight =  -std::numeric_limits<double>::max();
  unsigned int n_elem = PF_cloud.size();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:max_weight)
#endif
  for(unsigned int i = 0; i < n_elem; ++i){ // loop over cloud elements
    particle &p = PF_cloud[i];
    std::unique_ptr<dist_comb> &dc = ans.dists[i];

    double log_prob_y_given_state = y_dist->log_dens(dc->get_mean()),
      log_prop_transition = is_forward ?
        dens_calc.log_prob_state_given_parent(
          dc->get_mean(), p.get_state()) :
        dens_calc.log_prob_state_given_child(
          dc->get_mean(), p.get_state()),
      log_prop_proposal = dc->log_density(dc->get_mean());

    p.log_resampling_weight =
      p.log_weight + log_prop_transition + log_prob_y_given_state
      - log_prop_proposal;

    if(!is_forward)
      p.log_resampling_weight += (
        prior->log_dens(dc->get_mean()) - prior_p1->log_dens(p.get_state()));

    max_weight = MAX(max_weight, p.log_resampling_weight);
  } // end loop over cloud elements

  auto norm_out =
    normalize_log_resampling_weight
    <true, true>
    (PF_cloud, max_weight);
  outcome = sample(data, norm_out.weights, norm_out.ESS, did_resample);

  return std::move(ans.dists);
}

template class AUX_resampler_normal_approx_w_particles<false>;
template class AUX_resampler_normal_approx_w_particles<true>;
