#ifndef RESAMPLERS
#define RESAMPLERS

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"

/*
 Each "resampler" has a static function which:
    1) computes log re-sampling weights assuming that weights are computed
    2) samples according to re-sampling weights
    3) returns re-sampled indices if effective sample size is low. Otherwise
       return an index for each element in the cloud
*/

/*
 Non-auxiliary particle filter where we just to the weights as they are. I.e.
 set:
  beta_j = w_{j - 1}
*/
template<typename densities>
class None_AUX_resampler {
public:
  inline static arma::uvec resampler(PF_data &data, cloud &PF_cloud, unsigned int t){

    /* Compute effective sample size (ESS) */
    arma::vec weights(data.N_AUX);
    double ESS;
    auto w = weights.begin();
    for(auto it = PF_cloud.begin(); it != PF_cloud.end(); ++it, ++w){
      /* No need to update weights */
      it->log_resampling_weight = it->log_weight;

      *w = exp(it->log_resampling_weight);
      ESS += *w * *w;
    }
    ESS = 1/ ESS;

    if(ESS < data.forward_backward_ESS_threshold){
      return(sample_indices(weights));
    }

    return(arma::linspace<arma::uvec>(0, data.N_AUX - 1, data.N_AUX));
  }
};

#endif
