#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

#define MAX(a,b) (((a)>(b))?(a):(b))

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"
#include "PF_utils.h"

/*
  Each "importance_sampler" has two static functions:
    sample:
      Returns a new cloud of particles with sampled states. The states are
      sampled acording to the specific importance density. The function also
      sets log_importance_dens on each particle
    sample_smooth:
      Returns sample given past and next state
    log_importance_dens_smooth:
      Returns the log importance density for state smoothing
    sample_first_state_n_set_weights:
      Returns a particle cloud for time zero or d + 1 with weights set
*/

/* base class importance samplers */

template<typename densities, bool is_forward>
class importance_dens_base {
public:
  static cloud sample_first_state_n_set_weights(const PF_data &data){
    cloud ans;
    ans.reserve(data.N_first);
    arma::mat m1;
    const arma::mat *Q_chol;
    const arma::vec *mean = &data.a_0;

    if(is_forward){
      Q_chol = &data.Q_0.chol;

    } else {
      m1 = arma::chol(data.Q.mat * data.d + data.Q_0.mat);
      Q_chol = &m1;
    }

    if(data.debug > 1){
      data.log(2) << "Sampling "
                  << (is_forward ? "first" : "state d + 1")
                  << " with chol(covariance matrix):";
      data.log(2) << *Q_chol;
      data.log(2) << "and mean:";
      data.log(2) << mean->t();
    }

    double log_weight = log(1. / data.N_first);
    for(arma::uword i = 0; i < data.N_first; ++i){
      arma::vec new_state = mvrnorm(*mean, *Q_chol);
      ans.New_particle(new_state, nullptr);
      ans[i].log_weight = log_weight;
    }

    return(ans);
  }
};

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

 See:
  Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class importance_dens_no_y_dependence :
  public importance_dens_base<densities, is_forward>{
  static double log_importance_dens(const PF_data &data, const particle &p, int t){
    /* independent of is_forward for first order random walk */
    return dmvnrm_log(p.state, p.parent->state, data.Q_proposal.chol_inv);
  }

  static double log_importance_dens_smooth(
      const PF_data &data, const particle &p, int t){
    arma::vec mean = p.parent->state + p.child->state;
    mean *= .5;

    return dmvnrm_log(p.state, mean, data.Q_proposal_smooth.chol_inv);
  }

public:
  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t, nothing unused){
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by";
      data.log(3) << data.Q_proposal.chol;
    }

    auto it = resample_idx.begin();
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      arma::vec new_state = mvrnorm(cl[*it].state, data.Q_proposal.chol);
      ans.New_particle(new_state, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = log_importance_dens(data, p, t);
    }

    return ans;
  }

  static cloud sample_smooth(
    const PF_data &data,
    cloud &fw_cloud, const arma::uvec &fw_idx,
    cloud &bw_cloud, const arma::uvec &bw_idx, const unsigned int t){
    cloud ans;
    ans.reserve(data.N_smooth);

    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by";
      data.log(3) << data.Q_proposal_smooth.chol;
    }

    auto it_fw = fw_idx.begin();
    auto it_bw = bw_idx.begin();
    for(arma::uword i = 0; i < data.N_smooth; ++i, ++it_fw, ++it_bw){
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mean = fw_p.state + bw_p.state;
      mean *= .5;

      arma::vec new_state = mvrnorm(
        mean, data.Q_proposal_smooth.chol);
      ans.New_particle(new_state, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = log_importance_dens_smooth(data, p, t);
    }

    return ans;
  }
};

/*
  Sampler which makes a normal approximation for the observed outcome. See the
  second example on page 462-463 of:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class importance_dens_normal_approx  :
  public importance_dens_base<densities, is_forward> {
  inline static void debug_msg_before_sampling(
      const PF_data &data,
      const input_for_normal_approximation &inter_output){
    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Sigma) given by";
      data.log(3) << inter_output.Sigma_chol;
      data.log(3) << "The mean before accounting for the parent (and child) particle is:";
      data.log(3) << solve_w_precomputed_chol(inter_output.Sigma_inv_chol, inter_output.mu).t();
    }
  }

  inline static void debug_msg_while_sampling(
      const PF_data &data, const particle &p, const arma::vec &mu){
    if(data.debug > 4){
      data.log(5) << "Sampled particle:";
      data.log(5) << p.state.t();
      data.log(5) << "from normal distribution with mean:";
      data.log(5) << mu.t();
      data.log(5) << "The parent had state:";
      data.log(5) << p.parent->state.t();

      if(p.child){
        data.log(5) << "and the child had state";
        data.log(5) << p.child->state.t();
      }
    }
  }

public:
  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t, nothing unused){
    /* Find weighted mean estimate */
    arma::vec alpha_bar = cl.get_weigthed_mean();

    /* compute parts of the terms for the mean and covariance */
    auto &Q = data.Q_proposal;
    auto inter_output = compute_mu_n_Sigma_from_normal_approximation<densities>(
      data, t, Q, alpha_bar, cl);

    return(sample(data, cl, resample_idx, t, inter_output));
  }

  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t,
      input_for_normal_approximation &inter_output){
    /* Sample */
    debug_msg_before_sampling(data, inter_output);

    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    auto it = resample_idx.begin();
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      arma::vec &mu_j = inter_output.mu_js[*it];

      arma::vec new_state = mvrnorm(mu_j, inter_output.Sigma_chol);
      ans.New_particle(new_state, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.state, mu_j, inter_output.sigma_chol_inv);

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }

  static cloud sample_smooth(
      const PF_data &data,
      cloud &fw_cloud, const arma::uvec &fw_idx,
      cloud &bw_cloud, const arma::uvec &bw_idx, const unsigned int t){
    /* Find weighted mean estimate */
    arma::vec alpha_bar = fw_cloud.get_weigthed_mean() + bw_cloud.get_weigthed_mean();
    alpha_bar *= .5;

    /* compute parts of the terms for the mean and covariance */
    auto &Q = data.Q_proposal_smooth;
    auto inter_output = compute_mu_n_Sigma_from_normal_approximation<densities>(
      data, t, Q, alpha_bar);

    /* Sample */
    debug_msg_before_sampling(data, inter_output);

    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    auto it_fw = fw_idx.begin();
    auto it_bw = bw_idx.begin();
    for(arma::uword i = 0; i < data.N_smooth; ++i, ++it_fw, ++it_bw){
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mu_j = fw_p.state + bw_p.state;
      mu_j *= .5;
      mu_j = solve_w_precomputed_chol(Q.chol, mu_j) + inter_output.mu;
      mu_j = solve_w_precomputed_chol(inter_output.Sigma_inv_chol, mu_j);

      arma::vec new_state = mvrnorm(mu_j, inter_output.Sigma_chol);
      ans.New_particle(new_state, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.state, mu_j, inter_output.sigma_chol_inv);

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }
};



#undef MAX
#endif
