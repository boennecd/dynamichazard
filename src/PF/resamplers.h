#ifndef RESAMPLERS
#define RESAMPLERS

#include "PF_utils.h"
#include "../sample_funcs.h"
#include "../utils.h"

/* base class for common functions*/

template <arma::uvec (*sample_func)(const arma::uword, arma::vec&)>
class resampler_base {
protected:
  static arma::uvec sample(const PF_data&, arma::vec&, const double, bool&);
};

/*
 Each "resampler" has a static function which:
    1) computes log re-sampling weights assuming that weights are computed
    2) samples according to re-sampling weights
    3) returns re-sampled indices if effective sample size is low through the
       `outcome` argument. The actual returned object contains intermediate
       object to avoid recomputation in the importance sampler. If the
       effective sample size i large than the `outcome` argument is  an index
       for each element in the cloud. The boolean argument is used to indicate
       if sampling is performed.
*/

/*
 Non-auxiliary particle filter where we just to the weights as they are. E.g.
 set
    beta_j = w_{j,t - 1}
 in the forward particle filter.
*/

template<bool is_forward>
class None_AUX_resampler : private resampler_base<systematic_resampling> {
public:
  static nothing resampler(
      pf_dens&, const PF_data&, cloud&, std::shared_ptr<PF_cdist>,
      arma::uword, arma::uvec &, bool&);
};

/*
  Auxiliary particle filter with weights as in the end of page 462 of:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing
    algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<bool is_forward>
class AUX_resampler_normal_approx_w_cloud_mean :
  private resampler_base<systematic_resampling> {
public:
  static std::vector<std::unique_ptr<dist_comb>> resampler(
      pf_dens&, const PF_data&, cloud&, std::shared_ptr<PF_cdist>,
      arma::uword, arma::uvec &, bool&);
};

/*
  Auxiliary particle filter with weights as in the end of page 462 of the
  following paper with Taylor expansion around the parent particle:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing
    algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<bool is_forward>
class AUX_resampler_normal_approx_w_particles :
  private resampler_base<systematic_resampling> {
public:
  static std::vector<std::unique_ptr<dist_comb>> resampler(
      pf_dens&, const PF_data&, cloud&, std::shared_ptr<PF_cdist>,
      arma::uword, arma::uvec &, bool&);
};

#endif
