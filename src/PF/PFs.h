#ifndef AUX_PF_H
#define AUX_PF_H

#define BOOT_FILTER "bootstrap_filter"
#define PF_APPROX_CLOUD_MEAN "PF_normal_approx_w_cloud_mean"
#define AUX_APPROX_CLOUD_MEAN "AUX_normal_approx_w_cloud_mean"
#define PF_APPROX_PARTICLE "PF_normal_approx_w_particles"
#define AUX_APPROX_PARTICLE "AUX_normal_approx_w_particles"

#include "PF_utils.h"

class PF_base {
protected:
  static void debug_msg_after_weighting(
      const PF_data&, cloud&, const bool have_resampled = false,
      const unsigned int max_size = 0L);
};

/*
  Auxiliary particle filter as described in:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing
    algorithm with linear computational cost. Biometrika, 97(2), 447-464.

  The forward filter returns a particle cloud for time 0, 1, ..., d
  The backward filter returns a particle cloud for time d + 1, d, ..., 1
*/

template<
    template <bool> class T_resampler,
    template <bool> class T_importance_dens,
    bool is_forward>
class AUX_PF : private PF_base {
  using resampler = T_resampler<is_forward>;
  using importance_dens = T_importance_dens<is_forward>;

public:
  static std::vector<cloud>
  compute(const PF_data&, pf_dens&);
};

/*
  O(N) smoother from:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing
    algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<
  template <bool> class T_resampler,
  template <bool> class T_importance_dens>
class PF_smoother_Fearnhead_O_N : private PF_base {
  using uword = arma::uword;
  static arma::uvec sample_idx(const PF_data&, cloud&);

public:
  static smoother_output
  compute(const PF_data&, pf_dens&);

  using forward_filter =
    AUX_PF<T_resampler, T_importance_dens, true>;
  using backward_filter =
    AUX_PF<T_resampler, T_importance_dens, false>;
  using importance_dens = T_importance_dens<false /* arg should not matter*/>;
};

/*
  O(N^2) smoother from:
    Briers, M., Doucet, A., & Maskell, S. (2010). Smoothing algorithms for
    state–space models. Annals of the Institute of Statistical Mathematics,
    62(1), 61-89.
*/

template<
  template <bool> class T_resampler,
  template <bool> class T_importance_dens>
class PF_smoother_Brier_O_N_square : private PF_base {
  using uword = arma::uword;

public:
  static smoother_output
  compute(const PF_data&, pf_dens&);

  using forward_filter =
    AUX_PF<T_resampler, T_importance_dens, true>;
  using backward_filter =
    AUX_PF<T_resampler, T_importance_dens, false>;
};

#endif
