#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"

/*
  Each "importance_sampler" has two static functions:
    sample:
      Returns a new cloud of particles with sampled states. The states are
      sampled acording to the specific importance density
    log_importance_dens:
      Returns the log importance density
    sample_smooth:
      Returns sample given past and next state
    log_importance_dens_smooth:
      Returns the log importance density for state smoothing
    sample_state_zero_n_set_weights:
      Returns a particle cloud for time zero with weights set
    sample_d_plus_1_state_n_set_weights:
      Returns a particle cloud for time d + 1 with weights set
*/

/*
 Importance sampler with importance density that does not depend on the
 outcome for the first order random walk. That is:
  q(alpah_t | alpha_{t - 1}^{j}, y_t) =
    N(alpha | alpha_{t - j}^{j}, Q)
  tilde{q}(alpah_t | alpha_{t + 1}^{k}, y_t) =
    N(alpha | alpha_{t + 1}^{j}, Q)
  tilde{q}(alpah_t | alpha_{t + 1}^{j}, alpha_{t + 1}^{k}, y_t) =
    N(alpha | (1/2)Q(Q^-1 alpha_{t + 1}^{j} + alpha_{t + 1}^{k}), (1/2)Q) =
    N(alpha | (1/2)alpha_{t + 1}^{j} + (1/2)alpha_{t + 1}^{k}), (1/2)Q)

 See page 461-462 of:
  Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/
template<typename densities, bool is_forward>
class importance_dens_no_y_dependence {
public:
  static cloud sample(const PF_data &data, cloud &cl, const arma::uvec &resample_idx, const unsigned int t){
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    auto it = resample_idx.begin();
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      arma::vec new_state = mvrnorm(cl[*it].state, data.Q_chol);
      ans.New_particle(new_state, &cl[*it]);
    }

    return ans;
  }

  static double log_importance_dens(const PF_data &data, const particle &p, int t){
    /* independent of is_forward for first order random walk */
    return dmvnrm_log(p.state, p.parent->state, data.Q_chol_inv);
  }

  static cloud sample_smooth(
    const PF_data &data,
    cloud &fw_cloud, const arma::uvec &fw_idx,
    cloud &bw_cloud, const arma::uvec &bw_idx, const unsigned int t){
    cloud ans;
    ans.reserve(data.N_smooth);

    auto it_fw = fw_idx.begin();
    auto it_bw = fw_idx.begin();
    for(arma::uword i = 0; i < data.N_smooth; ++i, ++it_fw, ++it_bw){
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec new_state = mvrnorm(
        1/2 * (fw_p.state + bw_p.state), data.Q_half_chol);
      ans.New_particle(new_state, &fw_p, &bw_p);
    }

    return ans;
  }

  static double log_importance_dens_smooth(
      const PF_data &data, const particle &p, int t){
    return dmvnrm_log(p.state, 1/2 * (p.parent->state + p.child->state), data.Q_half_chol_inv);
  }

  static cloud sample_state_zero_n_set_weights(const PF_data &data){
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    double log_weight = log(data.N_fw_n_bw);
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      // Do not sample but set to a_0 with equal weight on each
      ans.New_particle(data.a_0, nullptr);
      ans[i].log_weight = log_weight;
    }

    return(ans);
  }

  static cloud sample_d_plus_1_state_n_set_weights(const PF_data &data){
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    const arma::mat Q_d_chol = data.Q_chol * sqrt(data.d + 1);
    const arma::mat Q_chol_inv = arma::inv(arma::trimatu(Q_d_chol));
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      arma::vec new_state = mvrnorm(data.a_0, Q_d_chol);
      ans.New_particle(new_state, nullptr);
      ans[i].log_weight = dmvnrm_log(new_state, data.a_0, Q_chol_inv);
    }

    return(ans);
  }

};

#endif
