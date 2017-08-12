#ifndef RESAMPLERS
#define RESAMPLERS

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"
#include "PF_utils.h"


#define MAX(a,b) (((a)>(b))?(a):(b))

/* base class for common functions*/

template <arma::uvec (*sample_func)(const arma::uword, arma::vec&)>
class resampler_base {
protected:
  inline static arma::uvec sample(
      const PF_data &data, arma::vec &probs, const double ESS){
    if(probs.n_elem != data.N_fw_n_bw){
      if(data.debug > 1){
        data.log(2) << "Subsampling " << probs.n_elem << " to get "
                    << data.N_fw_n_bw <<  " using re-sampling weights";
      }

      return(sample_func(data.N_fw_n_bw, probs));
    }

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

      return(sample_func(data.N_fw_n_bw, probs));
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
class None_AUX_resampler : private resampler_base<systematic_resampling> {
public:
  inline static arma::uvec resampler(const PF_data &data, cloud &PF_cloud, unsigned int t){
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

    return(sample(data, weights, ESS));
  }
};

/*
  Auxiliary particle filter where we use the weights times the next period
  outcome where we just use the previous state. I.e.:
    beta_j = w_{j - 1} * P(y_t|alpha_{t - 1})
  or
    beta_j = w_{j - 1} * P(y_t|alpha_{t + 1})
  where t depends on the direction
*/

template<typename densities, bool is_forward, bool scale_w_proposal>
class AUX_resampler_tp : private resampler_base<systematic_resampling> {
public:
  inline static arma::uvec resampler(const PF_data &data, cloud &PF_cloud, unsigned int t){
    /* Compute effective sample size (ESS) */
    arma::vec weights(PF_cloud.size());

    /* compute non-normalized log weights */
    bool parent_has_log_importance_dens =
      !((t == 1 && is_forward) || (t == (unsigned)data.d && !is_forward));
    double max_weight =  -std::numeric_limits<double>::max();
    for(auto it = PF_cloud.begin(); it != PF_cloud.end(); ++it){
      double log_prob_y_given_state = densities::log_prob_y_given_state(
        data, *it /* will either be state from t - 1 or t + 1*/, t);

      it->log_resampling_weight = it->log_weight + log_prob_y_given_state;
      if(scale_w_proposal){
        if(parent_has_log_importance_dens){
          it->log_resampling_weight -= it->log_importance_dens;
        }
      }

      max_weight = MAX(it->log_resampling_weight, max_weight);
    }

    auto norm_out = normalize_log_resampling_weight<true, true>(PF_cloud, max_weight);

    return(sample(data, norm_out.weights, norm_out.ESS));
  }
};

template<typename densities, bool is_forward>
using AUX_resampler_no_proposal_scaling = AUX_resampler_tp<densities, is_forward, false>;

/*
 Auxiliary particle filter where we use the weights times the next period
 outcome where we just use the previous state and we devide by the propsal
 likelihood. I.e.:
  beta_j = w_{j - 1} * P(y_t|alpha_{t - 1}) / q(alpah_{j - 1} | alpha_{t - 2}, y_{t - 1})
 or
  beta_j = w_{j - 1} * P(y_t|alpha_{t + 1}) / q(alpah_{j + 1} | alpha_{t + 2}, y_{t + 1})
 where t depends on the direction
*/

template<typename densities, bool is_forward>
using AUX_resampler_with_proposal_scaling = AUX_resampler_tp<densities, is_forward, true>;

#undef MAX
#endif
