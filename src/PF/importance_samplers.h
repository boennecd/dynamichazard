#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"

/*
 Importance sampler with importance density that does not depend on the
 outcome. That is:
 q(alpah_t | ... ) ~ N(alpah_t | alpha_{parent}^{(j)}, Q)
 */
template<typename densities>
class importance_dens_n_y_dependence {
public:
  static cloud sample(const PF_data &data, const cloud &cl, arma::uvec resample_idx, const unsigned int t){
    cloud ans;

    for(auto it = resample_idx.begin(); it != resample_idx.end(); ++it){
      arma::vec new_state = mvrnorm(cl[*it].state, data.Q_chol);
      ans.New_particle(new_state, cl[*it]);
    }

    return ans;
  }
};

#endif
