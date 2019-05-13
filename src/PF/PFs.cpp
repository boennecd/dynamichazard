#include "PFs.h"
#include "importance_samplers.h"
#include "resamplers.h"
#include "../sample_funcs.h"
#include "../utils.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void PF_base::debug_msg_after_weighting(
    const PF_data &data, cloud &cl, const bool have_resampled,
    const unsigned int max_size){
  if(data.debug > 1){
    double max_w =  -std::numeric_limits<double>::max();
    double min_w = std::numeric_limits<double>::max();
    double ESS = 0.;
    for(auto it = cl.begin(); it != cl.end(); ++it){
      max_w = MAX(max_w, it->log_weight);
      min_w = MIN(min_w, it->log_weight);

      double w = exp(it->log_weight);
      ESS += w * w;
    }
    ESS = 1 / ESS;

    if(have_resampled)
      data.log(2) << "Sub-sampled cloud. There are " << cl.size()
                  << " unique particles where up to " << max_size
                  << " is possible. ";

    data.log(2) << "(min, max) log weights are: ("
                << min_w  << ", " << max_w  <<  "). "
                << "ESS (before re-weighting) is: " << ESS;
  }
}



template<
  template <bool> class T_resampler,
  template <bool> class T_importance_dens,
  bool is_forward>
std::vector<cloud>
AUX_PF<T_resampler, T_importance_dens, is_forward>::
  compute(const PF_data &data, pf_dens &dens_calc)
{
  std::vector<cloud> clouds;
  std::string direction_str = (is_forward) ? "forward" : "backward";

  if(data.debug > 0)
    data.log(1) << "Running " << direction_str << " filter"
                << "\nSampling first particle at time "
                << static_cast<std::string>(is_forward ? "0" : "d + 1");
  clouds.push_back(
    importance_dens::sample_first_state_n_set_weights(dens_calc, data));

  int t = is_forward ? 1 : data.d;
  for(int iter = 1; iter <= data.d; ++iter){
    if((iter + 1) % 3 == 0)
      Rcpp::checkUserInterrupt();

    std::shared_ptr<PF_cdist> y_dist = dens_calc.get_y_dist(t),
      prior, prior_p1;

    if(!is_forward){
      prior = dens_calc.get_prior(t);
      prior_p1 = dens_calc.get_prior(t + 1L);
    }

    /* re-sample indicies */
    if(data.debug > 0)
      data.log(1) << "Starting iteration " << t << ". Re-sampling weights";
    arma::uvec resample_idx;
    bool did_resample;
    auto additional_resampler_out = resampler::resampler(
      dens_calc, data, clouds.back(), y_dist, t, resample_idx, did_resample);

    if(data.debug > 0){
      if(did_resample){
        data.log(1) << "Did resample";
      } else
        data.log(1) << "Did not re-sample";
    }

    /* sample new cloud */
    if(data.debug > 0)
      data.log(1) << "Sampling states";
    cloud new_cloud = importance_dens::sample(
      y_dist, dens_calc, data, clouds.back(), resample_idx, t,
      additional_resampler_out);

    /* update weights */
    if(data.debug > 0)
      data.log(1) << "Updating weights";
    {
        double max_weight =  -std::numeric_limits<double>::max();
        arma::uvec r_set = get_risk_set(data.risk_sets, t);
        unsigned int n_elem = new_cloud.size();
        double log_N = std::log(n_elem);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:max_weight)
#endif
        for(unsigned int i = 0; i < n_elem; ++i){ // loop over new particles
          auto it = new_cloud.begin() + i;
          double log_prob_y_given_state =
            y_dist->log_dens(it->get_state());
          double log_prob_state_given_other =
            is_forward ?
            dens_calc.log_prob_state_given_parent(
              it->get_state(), it->parent->get_state()) :
            dens_calc.log_prob_state_given_child(
              it->get_state(), it->parent->get_state());

          it->log_likelihood_term = it->log_weight =
            /* nominator */
            log_prob_y_given_state + log_prob_state_given_other
            /* denoninator */
            - it->log_importance_dens;

            if(did_resample){
              it->log_weight +=
                it->parent->log_weight - it->parent->log_resampling_weight;
              /* See
              *   Doucet, A., & Johansen, A. M. (2009). A tutorial on particle
              *   filtering and smoothing: Fifteen years later. Handbook of
              * nonlinear filtering, 12(656-704), 3.
              * Page 11, 15, 21, and 26.
              */
              it->log_likelihood_term +=
                it->parent->log_weight - it->parent->log_resampling_weight -
                log_N;

            } else {
              it->log_weight += it->parent->log_weight;
              it->log_likelihood_term += it->parent->log_weight;

            }

            if(!is_forward){
              it->log_weight +=
                prior->log_dens(it->get_state()) -
                prior_p1->log_dens(it->parent->get_state());
            }

            max_weight = MAX(max_weight, it->log_weight);

        } // end loop over new particle

        normalize_log_weights<false, true>(new_cloud, max_weight);
    }

    debug_msg_after_weighting(data, new_cloud);

    /* Add cloud  */
    clouds.push_back(std::move(new_cloud));

    if(is_forward)
      ++t;
    else
      --t;
  }

  return(clouds);
}



template<
  template <bool> class T_resampler,
  template <bool> class T_importance_dens>
arma::uvec
PF_smoother_Fearnhead_O_N<T_resampler, T_importance_dens>::
  sample_idx(const PF_data &data, cloud &cl){
  auto size = cl.size();
  arma::vec probs(size);

  auto pr = probs.begin();
  auto part = cl.begin();
  for(uword j = 0; j < size; ++j, ++pr, ++part)
    *pr = exp(part->log_resampling_weight);

  return systematic_resampling(data.N_smooth, probs);
}



template<
  template <bool> class T_resampler,
  template <bool> class T_importance_dens>
smoother_output
PF_smoother_Fearnhead_O_N<T_resampler, T_importance_dens>::
  compute(const PF_data &data, pf_dens &dens_calc){
    smoother_output result;
    std::vector<cloud> &forward_clouds  = result.forward_clouds;
    std::vector<cloud> &backward_clouds = result.backward_clouds;
    std::vector<cloud> &smoothed_clouds = result.smoothed_clouds;

    forward_clouds  =  forward_filter::compute(data, dens_calc);
    backward_clouds = backward_filter::compute(data, dens_calc);

    if(data.debug > 0)
      data.log(1) << "Finished finding forward and backward clouds. Started smoothing";

    auto fw_cloud = forward_clouds.begin(); // first index is time 0
    auto bw_cloud = backward_clouds.rbegin();
    ++bw_cloud; // first index is time 1 -- we need index 2 to start with

    for(int t = 1; t <= data.d /* note the leq */;
    ++t, ++fw_cloud, ++bw_cloud){
    std::shared_ptr<PF_cdist> y_dist = dens_calc.get_y_dist(t),
      prior_p1 = dens_calc.get_prior(t + 1L);

      // TODO: maybe use backward cloud at time t == 1
      if(t == data.d){ // we can use the forward particle cloud
        ++fw_cloud; // need time d cloud
        cloud last_cloud = *fw_cloud; // copy

        debug_msg_after_weighting(data, last_cloud);

        smoothed_clouds.push_back(std::move(last_cloud));

        continue;
      }

      /* re-sample */
      if(data.debug > 0)
        data.log(1) << "Started smoothing at time " << t
                    << "\nRe-sampling indices of previous and next state";

      arma::uvec fw_idx = sample_idx(data, *fw_cloud); // sample forward particles
      arma::uvec bw_idx = sample_idx(data, *bw_cloud); // sample backward particle

      /* sample states */
      if(data.debug > 0)
        data.log(1) << "Sampling states of previous and next state";

      cloud new_cloud = importance_dens::sample_smooth(
        y_dist, dens_calc, data, *fw_cloud, fw_idx, *bw_cloud, bw_idx, t);

      /* update weight */
      if(data.debug > 0)
        data.log(1) << "Weighting particles";
      {
        double max_weight = -std::numeric_limits<double>::max();
        arma::uvec r_set = get_risk_set(data.risk_sets, t);
        unsigned int n_elem = new_cloud.size();

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:max_weight)
#endif
        for(unsigned int i = 0; i < n_elem; ++i){ // loop over smooth clouds
          auto it = new_cloud.begin() + i;
          double log_prob_y_given_state = y_dist->log_dens(it->get_state()),
            log_prob_state_given_previous =
              dens_calc.log_prob_state_given_parent(
                it->get_state(), it->parent->get_state()),
                log_prob_next_given_state =
                  dens_calc.log_prob_state_given_child(
                    it->get_state(), it->child->get_state()),
                    log_importance_dens = it->log_importance_dens,
                    log_artificial_prior = // TODO: have already been computed
                      prior_p1->log_dens(it->child->get_state());

          it->log_likelihood_term = it->log_weight =
            /* nominator */
            (log_prob_y_given_state + log_prob_state_given_previous + log_prob_next_given_state +
            it->parent->log_weight + it->child->log_weight)
            /* denoninator */
            - (log_importance_dens + it->parent->log_resampling_weight +
            it->child->log_resampling_weight + log_artificial_prior);

          max_weight = MAX(max_weight, it->log_weight);
        } // end over smooth clouds

        normalize_log_weights<false, true>(new_cloud, max_weight);
      }

      debug_msg_after_weighting(data, new_cloud);

      if(data.N_smooth_final < data.N_smooth){
        new_cloud = re_sample_cloud(data.N_smooth_final, new_cloud);
        debug_msg_after_weighting(data, new_cloud, true, data.N_smooth_final);
      }

      /* Add cloud  */
      smoothed_clouds.push_back(std::move(new_cloud));
    }

    return result;
  }

template class PF_smoother_Fearnhead_O_N<
  None_AUX_resampler,
  importance_dens_no_y_dependence>;
template class PF_smoother_Fearnhead_O_N<
  None_AUX_resampler,
  importance_dens_normal_approx_w_cloud_mean>;
template class PF_smoother_Fearnhead_O_N<
  AUX_resampler_normal_approx_w_cloud_mean,
  importance_dens_normal_approx_w_cloud_mean>;
template class PF_smoother_Fearnhead_O_N<
  None_AUX_resampler,
  importance_dens_normal_approx_w_particles>;
template class PF_smoother_Fearnhead_O_N<
  AUX_resampler_normal_approx_w_particles,
  importance_dens_normal_approx_w_particles>;



template<
  template <bool> class T_resampler,
  template <bool> class T_importance_dens>
smoother_output
PF_smoother_Brier_O_N_square<T_resampler, T_importance_dens>::
  compute(const PF_data &data, pf_dens &dens_calc){
    smoother_output result;
    std::vector<cloud> &forward_clouds  = result.forward_clouds;
    std::vector<cloud> &backward_clouds = result.backward_clouds;
    std::vector<cloud> &smoothed_clouds = result.smoothed_clouds;
    std::shared_ptr<smoother_output::trans_like_obj> trans_ptr =
      result.get_transition_likelihoods(false);

    smoother_output::trans_like_obj &transition_likelihoods = *trans_ptr;

    forward_clouds  = forward_filter::compute(data, dens_calc);
    backward_clouds = backward_filter::compute(data, dens_calc);

    if(data.debug > 0)
      data.log(1) << "Finished finding forward and backward clouds. Started smoothing";

    auto fw_cloud = forward_clouds.begin();
    //++fw_cloud; // first index is time 0
    auto bw_cloud = backward_clouds.rbegin();
    //++bw_cloud; // first index is time 1

    double max_weight = -std::numeric_limits<double>::max();
    for(int t = 1; t <= data.d /* note the leq */; ++t, ++fw_cloud, ++bw_cloud){
      // TODO: maybe use backward cloud at time t == 1
      if(t == data.d){ // we can use the forward particle cloud
        ++fw_cloud; // need time d cloud
        cloud last_cloud = *fw_cloud; // copy

        std::vector<smoother_output::particle_pairs> last_trans_like;
        last_trans_like.reserve(last_cloud.size());

        for(cloud::const_iterator pa = last_cloud.begin();
            pa != last_cloud.end(); ++pa){
          std::vector<smoother_output::pair> pairs(1);
          pairs[0].p = pa->parent;
          pairs[0].log_weight = 0;
          last_trans_like.emplace_back(
            &(*pa), pa->log_weight, std::move(pairs));
        }

        debug_msg_after_weighting(data, last_cloud);

        /* Add new elements */
        transition_likelihoods.push_back(std::move(last_trans_like));
        smoothed_clouds.push_back(std::move(last_cloud));

        continue;
      }

      std::shared_ptr<PF_cdist> prior = dens_calc.get_prior(t);

      if(data.debug > 0)
        data.log(1) << "Started smoothing at time " << t;

      unsigned int n_elem_fw = fw_cloud->size();
      unsigned int n_elem_bw = bw_cloud->size();

      std::vector<smoother_output::particle_pairs> new_trans_like(n_elem_bw);
      cloud new_cloud(n_elem_bw);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(max:max_weight)
#endif
      for(unsigned int i = 0; i < n_elem_bw; ++i){ // loop over bw particle
        particle &bw_particle = (*bw_cloud)[i];
        smoother_output::particle_pairs
          pf_pairs(&bw_particle, bw_particle.log_weight);
        pf_pairs.transition_pairs.reserve(n_elem_fw);

        double max_weight_inner = -std::numeric_limits<double>::max();
        for(auto fw_particle = fw_cloud->begin();
            fw_particle != fw_cloud->end();
            ++fw_particle){
          // compute un-normalized weight of pair
          double log_like_transition = dens_calc.log_prob_state_given_parent(
            bw_particle.get_state(), fw_particle->get_state());
          double this_term = fw_particle->log_weight + log_like_transition;

          // add pair information and update max log term seen so far
          pf_pairs.transition_pairs.emplace_back(&(*fw_particle), this_term);
          max_weight_inner = MAX(max_weight_inner, this_term);
        }

        //  compute log sum of logs and normalize weights of transition_pairs
        double this_log_weight;
        {
          auto tmp =
            normalize_log_weights
          <false, true>
          (pf_pairs.transition_pairs, max_weight_inner);
          this_log_weight = tmp.log_sum_logs;
        }

        // add pairs
        new_trans_like[i] = std::move(pf_pairs);

        double log_artificial_prior = // TODO: have already been computed
          prior->log_dens(bw_particle.get_state());

        // add likelihood from BW particle
        this_log_weight += bw_particle.log_weight - log_artificial_prior;

        // add particle to smooth cloud with weight
        particle &new_p = new_cloud.set_particle(i, bw_particle.get_state());
        new_p.log_likelihood_term = new_p.log_weight = this_log_weight;

        max_weight = MAX(max_weight, this_log_weight);
      } // end loop over bw particle

      // normalize smoothed weights
      normalize_log_weights<false, true>(new_cloud, max_weight);

      // Set particle on new_trans_like
      auto it_cl = new_cloud.begin();
      auto it_trans = new_trans_like.begin();
      for(unsigned int i = 0; i < n_elem_bw; ++i, ++it_cl, ++it_trans)
        it_trans->p = &(*it_cl);

      debug_msg_after_weighting(data, new_cloud);

      /* Add new elements */
      transition_likelihoods.push_back(std::move(new_trans_like));
      smoothed_clouds.push_back(std::move(new_cloud));
    }

    return result;
  }

template class PF_smoother_Brier_O_N_square<
    None_AUX_resampler,
    importance_dens_no_y_dependence>;
template class PF_smoother_Brier_O_N_square<
    None_AUX_resampler,
    importance_dens_normal_approx_w_cloud_mean>;
template class PF_smoother_Brier_O_N_square<
    AUX_resampler_normal_approx_w_cloud_mean,
    importance_dens_normal_approx_w_cloud_mean>;
template class PF_smoother_Brier_O_N_square<
    None_AUX_resampler,
    importance_dens_normal_approx_w_particles>;
template class PF_smoother_Brier_O_N_square<
    AUX_resampler_normal_approx_w_particles,
    importance_dens_normal_approx_w_particles>;
