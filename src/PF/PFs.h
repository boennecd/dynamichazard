#ifndef AUX_PF_H
#define AUX_PF_H

#include "PF_data.h"

template<
    template <typename> class T_resampler,
    template <typename> class T_importance_dens,
    class densities
>
class AUX_PF {
  using resampler = T_resampler<densities>;
  using importance_dens = T_importance_dens<densities>;

  std::vector<cloud> clouds;
  PF_data &data;

public:
  AUX_PF(PF_data &data): data(data) {}

  void compute(){
    /* TODO: what to do at time 0? */
    clouds.push_back(new cloud());
    cloud *current_cloud = clouds.back();

    for(unsigned int t = 0; t < data.d; ++t){
      /* re-sample indicies */
      arma::uvec resample_idx =
        resampler::resampler(data, clouds.back(), t);

      /* sample new cloud */
      cloud new_cloud =
        importance_dens::sample(data, clouds.back(), resample_idx, t);

      /* update weights */
      double max_Weigths =  std::numeric_limits<double>::min();
      for(auto it = new_cloud.begin(); it != new_cloud.end(); ++it){
        double log_prob_y_given_state =
          densities::log_prob_y_given_state(data, *it, t);
        double log_prob_state_given_parent =
          densities::log_prob_state_given_previous(data, *it, t);
        double log_importance_dens =
          importance_dens::log_importance_dens(data, *it, t);

        it->log_weight =
          /* nominator */
          (log_prob_y_given_state + log_prob_state_given_parent + it->parent->log_weight)
          /* denoninator */
          - (log_importance_dens + it->parent->log_resampling_weight);

        max_Weigths = std::max(it->log_weight, max_Weigths);
      }

      double norm_constant = 0;
      for(auto it = new_cloud.begin(); it != new_cloud.end(); ++it){
        /* back transform weights */
        it->log_weight = exp(it->log_weight - max_Weigths);
        norm_constant += it->log_weight;
      }

      for(auto it = new_cloud.begin(); it != new_cloud.end(); ++it){
        /* Re-scale and take log */
        it->log_weight = log(it->log_weight / norm_constant);
      }

      /* Add cloud  */
      clouds.push_back(std::move(new_cloud));
    }
  }
};

#endif
