#ifndef RESAMPLERS
#define RESAMPLERS

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"

class sample_wo_replacement{
protected:
  inline static arma::uvec sample(
      const PF_data &data, const arma::vec &probs, const double ESS){
    if(ESS < data.forward_backward_ESS_threshold){

      if(data.debug > 1){
        data.log(2) << "ESS is below threshold (" << ESS << " < "
                    << data.forward_backward_ESS_threshold << "). Re-sampling";
      }

      if(data.debug > 2){
        data.log(3) << "Re-sampling " << data.N_fw_n_bw << " indices "
                    << " from " <<  probs.n_elem << " elements "
                    << " with " <<  arma::max(probs) << " as the higest probability";
      }

      return(sample_indices(data.N_fw_n_bw, probs));
    }

    if(data.debug > 1){
      data.log(2) << "ESS is greater than threshold (" << ESS << " >= "
                  << data.forward_backward_ESS_threshold << "). No re-sampling needed";
    }

    return(arma::linspace<arma::uvec>(0, data.N_fw_n_bw - 1, data.N_fw_n_bw));
  }
};

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

template<typename densities, bool is_forward>
class None_AUX_resampler : private sample_wo_replacement {
public:
  inline static arma::uvec resampler(const PF_data &data, cloud &PF_cloud, unsigned int t){
    /* Compute effective sample size (ESS) */
    arma::vec weights(data.N_fw_n_bw);
    double ESS = 0;
    auto w = weights.begin();
    for(auto it = PF_cloud.begin(); it != PF_cloud.end(); ++it, ++w){
      /* No need to update weights */
      it->log_resampling_weight = it->log_weight;

      *w = exp(it->log_resampling_weight);
      ESS += *w * *w;
    }
    ESS = 1/ ESS;

    return(sample(data, weights, ESS));
  }
};

#endif
