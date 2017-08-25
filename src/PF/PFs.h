#ifndef AUX_PF_H
#define AUX_PF_H

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"
#include "PF_utils.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

class PF_base {
protected:
  static void debug_msg_after_weighting(const PF_data &data, cloud &cl){
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

      data.log(2) << "(min, max) log weights are: ("
                  << min_w  << ", " << max_w  <<  "). "
                  << "ESS (before re-weighting) is: " << ESS;
    }
  }
};

/*
  Auxiliary particle filter as described in:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.

  The forward filter returns a particle cloud for time 0, 1, ..., d
  The backward filter returns a particle cloud for time d + 1, d, ..., 1
*/

template<
    template <typename, bool> class T_resampler,
    template <typename, bool> class T_importance_dens,
    class densities,
    bool is_forward
>
class AUX_PF : private PF_base {
  using resampler = T_resampler<densities, is_forward>;
  using importance_dens = T_importance_dens<densities, is_forward>;

public:
  static std::vector<cloud> compute(const PF_data &data){
    std::vector<cloud> clouds;
    std::string direction_str = (is_forward) ? "forward" : "backward";

    /* TODO: what to do at time 0 or d + 1 */
    if(data.debug > 0)
      data.log(1) << "Running " << direction_str << " filter"
                  << "\nSampling first particle at time "
                  << static_cast<std::string>(is_forward ? "0" : "d + 1");
    clouds.push_back(importance_dens::sample_first_state_n_set_weights(data));

    int t = is_forward ? 1 : data.d;
    for(int iter = 1; iter <= data.d; ++iter){
      if((iter + 1) % 3 == 0)
        Rcpp::checkUserInterrupt();

      /* re-sample indicies */
      if(data.debug > 0)
        data.log(1) << "Starting iteration " << t << ". Re-sampling weights";
      arma::uvec resample_idx;
      bool did_resample;
      auto additional_resampler_out = resampler::resampler(data, clouds.back(), t, resample_idx, did_resample);

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
        data, clouds.back(), resample_idx, t, additional_resampler_out);

      /* update weights */
      if(data.debug > 0)
        data.log(1) << "Updating weights";
      {
          densities dens_calc = densities();
          double max_weight =  -std::numeric_limits<double>::max();
          for(auto it = new_cloud.begin(); it != new_cloud.end(); ++it){
            double log_prob_y_given_state =
              dens_calc.log_prob_y_given_state(data, *it, t);
            double log_prob_state_given_previous =
              is_forward ?
              dens_calc.log_prob_state_given_previous(
                data, it->state, it->parent->state, t) :
              /* Notice different order and t + 1 */
              dens_calc.log_prob_state_given_previous(
                data, it->parent->state, it->state, t + 1);

            it->log_unnormalized_weight = it->log_weight =
              /* nominator */
              log_prob_y_given_state + log_prob_state_given_previous
              /* denoninator */
              - it->log_importance_dens;

            if(did_resample){
              it->log_weight +=
                it->parent->log_weight - it->parent->log_resampling_weight;

            } else {
              it->log_weight += it->parent->log_weight;

            }

            if(!is_forward){
              it->log_weight += dens_calc.log_artificial_prior(data, *it, t);
              it->log_weight -= dens_calc.log_artificial_prior(data, *it->parent, t + 1);
            }

            max_weight = MAX(it->log_weight, max_weight);
          }

          normalize_log_weights<false, true>(new_cloud, max_weight);
      }

      debug_msg_after_weighting(data, new_cloud);

      /* Add cloud  */
      clouds.push_back(std::move(new_cloud));

      if(is_forward){
        ++t;
      } else
        --t;
    }

    return(clouds);
  }
};

/*
  O(N) smoother from:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<
  template <typename, bool> class T_resampler,
  template <typename, bool> class T_importance_dens,
  class densities
>
class PF_smoother : private PF_base {
  using uword = arma::uword;

  inline static arma::uvec sample_idx(const PF_data &data, cloud &cl){
    auto size = cl.size();
    arma::vec probs(size);

    auto pr = probs.begin();
    auto part = cl.begin();
    for(uword j = 0; j < size; ++j, ++pr, ++part)
      *pr = exp(part->log_resampling_weight);

    return systematic_resampling(data.N_smooth, probs);
  }

public:
  static smoother_output compute(const PF_data &data){
    smoother_output result;
    std::vector<cloud> &forward_clouds = result.forward_clouds;
    std::vector<cloud> &backward_clouds = result.backward_clouds;
    std::vector<cloud> &smoothed_clouds = result.smoothed_clouds;

    forward_clouds = forward_filter::compute(data);
    backward_clouds = backward_filter::compute(data);

    if(data.debug > 0)
      data.log(1) << "Finished finding forward and backward clouds. Started smoothing";

    auto fw_cloud = /* first index is time 0 */ forward_clouds.begin();
    auto bw_cloud = backward_clouds.rbegin();
    bw_cloud += 1; // first index is time 1 -- we need index 2 to start with

    for(int t = 1; t <= data.d /* note the leq */; ++t, ++fw_cloud, ++bw_cloud){
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
        data, *fw_cloud, fw_idx, *bw_cloud, bw_idx, t);

      /* update weight */
      if(data.debug > 0)
        data.log(1) << "Weighting particles";
      {
          densities dens_calc = densities();
          double max_weight = -std::numeric_limits<double>::max();
          for(auto it = new_cloud.begin(); it != new_cloud.end(); ++it){
            double log_prob_y_given_state =
              dens_calc.log_prob_y_given_state(data, *it, t);
            double log_prob_state_given_previous =
              dens_calc.log_prob_state_given_previous(
                data, it->state, it->parent->state, t);
            double log_prob_next_given_state =
              dens_calc.log_prob_state_given_previous(
                /* notice different order and t + 1*/
                data, it->child->state, it->state, t + 1);
            double log_importance_dens = it->log_importance_dens;
            double log_artificial_prior =
              dens_calc.log_artificial_prior(data, *it->child /* note child */, t + 1 /* note t + 1 */);

            it->log_unnormalized_weight = it->log_weight =
              /* nominator */
              (log_prob_y_given_state + log_prob_state_given_previous + log_prob_next_given_state +
                it->parent->log_weight + it->child->log_weight)
              /* denoninator */
              - (log_importance_dens + it->parent->log_resampling_weight +
                  it->child->log_resampling_weight + log_artificial_prior);


            max_weight = MAX(it->log_weight, max_weight);
          }

          normalize_log_weights<false, true>(new_cloud, max_weight);
      }

      debug_msg_after_weighting(data, new_cloud);

      /* Add cloud  */
      smoothed_clouds.push_back(std::move(new_cloud));
    }

    return result;
  }

  using forward_filter =
    AUX_PF<T_resampler, T_importance_dens, densities, true>;
  using backward_filter =
    AUX_PF<T_resampler, T_importance_dens, densities, false>;
  using importance_dens = T_importance_dens<densities, false /* arg should not matter*/>;
};

#undef MIN
#undef MAX

#endif
