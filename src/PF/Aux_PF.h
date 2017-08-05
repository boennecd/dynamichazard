#ifndef AUX_PF
#define AUX_PF

#include "../arma_n_rcpp.h"
#include "particles.h"

struct data_holder {
  const int d;
};

template<
    template <typename> class T_resampler,
    template <typename> class T_importance_dens,
    class obs_likelihood
>
class AUX_PF {
  using resampler = T_resampler<obs_likelihood>;
  using importance_dens = T_importance_dens<obs_likelihood>;

  std::vector<cloud> clouds;
  data_holder *const data;

public:
  AUX_PF(data_holder *data) data(data) {}

  void compute(){
    /* TODO: what to do at time 0? */
    clouds.add(new cloud());
    cloud *current_cloud = clouds.back();

    for(unsigned int t = 0; i < data.d; ++t){
      /* re-sample indicies */
      vector<int> resample_idx =
        resampler::resampler(data, clouds.back(), t);

      /* sample new cloud */
      cloud new_cloud =
        importance_dens::sample(data, clouds.back(), resample_idx, t);

      /* update weights */
      double max_Weigths =  std::numeric_limits<double>::min();
      for(auto it = new_cound.begin(); it != new_cound.end(); ++it){
        double log_prob_y_given_state =
          obs_likelihood::log_prob_y_given_state(data, *it, t);
        double log_prob_state_given_parent =
          obs_likelihood::log_prob_state_given_parent(data, *it, t);
        double log_importance_dens =
          importance_dens::log_importance_dens(data, *it, t);

        it->log_weight =
          /* nominator */
          (log_prob_y_given_state + log_prob_state_given_parent + it->parent->weight)
          /* denoninator */
          - (log_importance_dens + it->parent->resampling_weight);

        max_Weigths = std::max(it->log_weight, max_Weigths);
      }

      double norm_constant = 0;
      for(auto it = new_cound.begin(); it != new_cound.end(); ++it){
        /* back transform weights */
        it->log_weight = exp(it->log_weight - max_Weigths);
        norm_constant += it->log_weight;
      }

      for(auto it = new_cound.begin(); it != new_cound.end(); ++it){
        /* Re-scale and take log */
        it->log_weight = log(it->log_weight / norm_constant);
      }

      /* Add cloud  */
      clouds.push_back(std::move(new_cloud));
    }
  }
};

#endif
