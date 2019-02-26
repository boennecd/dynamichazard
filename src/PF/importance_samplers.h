#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

#include "PF_utils.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"

#define SAMPLE_SMOOTH_ARGS                                \
  std::shared_ptr<PF_cdist>, pf_dens&,                    \
  const PF_data&, cloud&, const arma::uvec&, cloud&,      \
  const arma::uvec&, const arma::uword

#define SAMPLE_COMMON_ARGS                                     \
  std::shared_ptr<PF_cdist>, pf_dens&, const PF_data&, cloud&, \
  const arma::uvec&, const arma::uword

/*
  Each "importance_sampler" has the following static functions:
    sample:
      Returns a new cloud of particles with sampled states. The states are
      sampled acording to the specific importance density. The function also
      sets log_importance_dens on each particle
    sample_smooth:
      Returns sample given past and next state
    sample_first_state_n_set_weights:
      Returns a particle cloud for time zero or d + 1 with weights set
*/

/* base class importance samplers */
template<bool is_forward>
class importance_dens_base {
public:
  static cloud sample_first_state_n_set_weights(pf_dens&, const PF_data&);
};

/*
 Bootstrap filter

 See:
  Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing
  algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<bool is_forward>
class importance_dens_no_y_dependence :
  public importance_dens_base<is_forward>{
public:
  static cloud sample(SAMPLE_COMMON_ARGS, nothing);
  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS);
};

/*
  Sampler which makes a normal approximation for the observed outcome made
  around the mean of the previous periods particle cloud. See the second
  example on page 462-463 of:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing
    algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<bool is_forward>
class importance_dens_normal_approx_w_cloud_mean  :
  public importance_dens_base<is_forward> {
  static void debug_msg_while_sampling(
      const PF_data&, const particle&, const arma::vec&);

public:
  static cloud sample(SAMPLE_COMMON_ARGS, nothing);
  static cloud sample
    (SAMPLE_COMMON_ARGS, std::vector<std::unique_ptr<dist_comb>>&);
  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS);
};

/*
 Same as above except that it is independent of the parent particle. See:
  Lin, Ming T., Junni L. Zhang, Qiansheng Cheng, and Rong Chen. "Independent
  particle filters." Journal of the American Statistical Association 100,
  no. 472 (2005): 1412-1421.
 Notice: the returned cloud will have null pointers to parents as they
         have no parents.
*/

template<bool is_forward>
class importance_dens_normal_approx_w_cloud_mean_independent  :
  public importance_dens_base<is_forward> {
  static void debug_msg_while_sampling(
      const PF_data&, const particle&, const arma::vec&);

public:
  static cloud sample(SAMPLE_COMMON_ARGS, nothing);
  static cloud sample
    (SAMPLE_COMMON_ARGS, std::unique_ptr<dist_comb>&);
};

/*
  Sampler which makes a normal approximation for the observed outcome made
  around the mean of the parent particle.
*/

template<bool is_forward>
class importance_dens_normal_approx_w_particles  :
  public importance_dens_base<is_forward> {
  static void debug_msg_while_sampling(
      const PF_data&, const particle&, const arma::vec&, const arma::mat&);

public:
  static cloud sample(SAMPLE_COMMON_ARGS, nothing);
  static cloud sample
    (SAMPLE_COMMON_ARGS, std::vector<std::unique_ptr<dist_comb>>&);
  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS);
};

#undef SAMPLE_SMOOTH_ARGS
#undef SAMPLE_COMMON_ARGS
#endif
