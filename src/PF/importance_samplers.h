#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

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
      m1 = arma::chol(densities::get_artificial_prior_covar(data, data.d));
      Q_chol = &m1;
    }

    if(data.debug > 1){
      data.log(2) << "Sampling "
                  << (is_forward ? "first" : "state d + 1")
                  << " with chol(covariance matrix):" << std::endl
                  << *Q_chol
                  << "and mean:" << std::endl
                  << mean->t();
    }

    double log_weight = log(1. / data.N_first);
    for(arma::uword i = 0; i < data.N_first; ++i){
      arma::vec new_state = mvrnorm(*mean, *Q_chol);
      ans.new_particle(new_state, nullptr);
      ans[i].log_weight = log_weight;
    }

    return(ans);
  }
};

/*
 Bootstrap filter

 See:
  Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class importance_dens_no_y_dependence :
  public importance_dens_base<densities, is_forward>{
  static double log_importance_dens_smooth(
      const PF_data &data, const particle &p, int t){
    arma::vec mean = p.parent->get_state() + p.child->get_state();
    mean *= .5;

    return dmvnrm_log(p.get_state(), mean, data.Q_proposal_smooth.chol_inv);
  }

public:
  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t, nothing unused){
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    const covarmat *Q_use;
    arma::vec *a_0_scaled;
    double bw_w1, bw_w2;
    if(is_forward){
      Q_use = &data.Q_proposal;

    } else {
      bw_w1 = (double)(t + 1) / (t + 2);
      bw_w2 =              1. / (t + 2);

      arma::vec &&tmp_vec = bw_w2 * data.a_0;
      a_0_scaled = new arma::vec(tmp_vec);
      arma::mat &&tmp = ((double)(t + 1) / (t + 2)) * data.Q.mat;
      Q_use = new covarmat(tmp);
    }

    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by" << std::endl
                  << data.Q_proposal.chol;
    }

    auto it = resample_idx.begin();
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      const arma::vec *mu;
      if(is_forward){
        mu = &cl[*it].get_state();

      } else {
        arma::vec &&tmp = bw_w1 * cl[*it].get_state() + *a_0_scaled;
        mu = new arma::vec(tmp);

      }

      arma::vec new_state = mvrnorm(*mu, Q_use->chol);
      ans.new_particle(new_state, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.get_state(), *mu, Q_use->chol_inv);

      if(!is_forward)
        delete mu;
    }

    if(!is_forward){
      delete Q_use;
      delete a_0_scaled;

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
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by" << std::endl
                  << data.Q_proposal_smooth.chol;
    }

    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mean = fw_p.get_state() + bw_p.get_state();
      mean *= .5;

      arma::vec new_state = mvrnorm(
        mean, data.Q_proposal_smooth.chol);
      ans.new_particle(std::move(new_state), &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = log_importance_dens_smooth(data, p, t);
    }

    return ans;
  }
};

/*
  Sampler which makes a normal approximation for the observed outcome made
  around the mean of the previous periods particle cloud. See the second
  example on page 462-463 of:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class importance_dens_normal_approx_w_cloud_mean  :
  public importance_dens_base<densities, is_forward> {
  inline static void debug_msg_before_sampling(
      const PF_data &data,
      const input_for_normal_apprx &inter_output){
    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Sigma) given by" << std::endl
                  << inter_output.Sigma_chol
                  << "The mean before accounting for the parent (and child) particle is:" << std::endl
                  << solve_w_precomputed_chol(inter_output.Sigma_inv_chol, inter_output.mu).t();
    }
  }

  inline static void debug_msg_while_sampling(
      const PF_data &data, const particle &p, const arma::vec &mu){
    if(data.debug > 4){
      auto log = data.log(5);
      log << "Sampled particle:" << std::endl
          << p.get_state().t()
          << "from normal distribution with mean:" << std::endl
          << mu.t()
          << "The parent had state:" << std::endl
          << p.parent->get_state().t();

      if(p.child){
        log << "and the child had state" << std::endl
            << p.child->get_state().t();
      }
    }
  }

public:
  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t, nothing unused){
    /* Find weighted mean estimate */
    arma::vec alpha_bar = cl.get_weigthed_mean();

    /* compute means and covariances */
    auto &Q = data.Q_proposal;
    auto inter_output = compute_mu_n_Sigma_from_normal_apprx_w_cloud_mean
      <densities, is_forward>
      (data, t, Q, alpha_bar, cl);

    return(sample(data, cl, resample_idx, t, inter_output));
  }

  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t,
      input_for_normal_apprx_w_cloud_mean &inter_output){
    /* Sample */
    debug_msg_before_sampling(data, inter_output);

    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      auto it = resample_idx.begin() + i;
      arma::vec &mu_j = inter_output.mu_js[*it];

      arma::vec new_state = mvrnorm(mu_j, inter_output.Sigma_chol);
      ans.new_particle(std::move(new_state), &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.get_state(), mu_j, inter_output.sigma_chol_inv);

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
    auto inter_output =
      compute_mu_n_Sigma_from_normal_apprx<densities, 2, true>(
        data, t, Q, alpha_bar);

    /* Sample */
    debug_msg_before_sampling(data, inter_output);

    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mu_j = fw_p.get_state() + bw_p.get_state();
      mu_j *= .5;
      mu_j = solve_w_precomputed_chol(Q.chol, mu_j) + inter_output.mu;
      mu_j = solve_w_precomputed_chol(inter_output.Sigma_inv_chol, mu_j);

      arma::vec new_state;
      new_state = mvrnorm(mu_j, inter_output.Sigma_chol);
      ans.new_particle(std::move(new_state), &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.get_state(), mu_j, inter_output.sigma_chol_inv);

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }
};

/*
  Sampler which makes a normal approximation for the observed outcome made
  around the mean of the parent particle.
*/

template<typename densities, bool is_forward>
class importance_dens_normal_approx_w_particles  :
  public importance_dens_base<densities, is_forward> {

  inline static void debug_msg_while_sampling(
      const PF_data &data, const particle &p, const arma::vec &mu, const arma::mat Sigma_chol){
    if(data.debug > 4){
      auto log = data.log(5);
      log << "Sampled particle:" <<  std::endl
          << p.get_state().t()
          << "from normal distribution with mean:"  <<  std::endl
          << mu.t()
          << "and chol(Sigma):"  <<  std::endl
          << Sigma_chol
          << "The parent had state:" <<  std::endl
          << p.parent->get_state().t();

      if(p.child){
        log << "and the child had state" <<  std::endl
            << p.child->get_state().t();
      }
    }
  }

public:
  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t, nothing unused){

    /* compute means and covariances */
    auto &Q = data.Q_proposal;
    auto inter_output =
      compute_mu_n_Sigma_from_normal_apprx_w_particles
      <densities, is_forward>
      (data, t, Q, cl);

    return(sample(data, cl, resample_idx, t, inter_output));
  }

  static cloud sample(
      const PF_data &data, cloud &cl, const arma::uvec &resample_idx,
      const unsigned int t,
      input_for_normal_apprx_w_particle_mean &inter_output){
    /* Sample */
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      auto it = resample_idx.begin() + i;
      auto &inter_o = inter_output[*it];

      arma::vec new_state;
      new_state = mvrnorm(inter_o.mu, inter_o.Sigma_chol);
      ans.new_particle(std::move(new_state), &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.get_state(), inter_o.mu, inter_o.sigma_chol_inv);

      debug_msg_while_sampling(data, p, inter_o.mu, inter_o.Sigma_chol);
    }

    return(ans);
  }

  static cloud sample_smooth(
      const PF_data &data,
      cloud &fw_cloud, const arma::uvec &fw_idx,
      cloud &bw_cloud, const arma::uvec &bw_idx, const unsigned int t){
    /* Compute means before accounting for outcomes */
    std::vector<arma::vec> mus(data.N_smooth);

    auto begin_fw = fw_idx.begin();
    auto begin_bw = bw_idx.begin();
    auto begin_mu = mus.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      arma::vec &mu_j = *(begin_mu + i);
      const particle &fw_p = fw_cloud[*(begin_fw + i)];
      const particle &bw_p = bw_cloud[*(begin_bw + i)];

      mu_j = fw_p.get_state() + bw_p.get_state();
      mu_j *= .5;
    }

    /* compute means and covariances */
    auto &Q = data.Q_proposal_smooth;
    auto inter_output =
      compute_mu_n_Sigma_from_normal_apprx_w_particles
      <densities, is_forward>
      (data, t, Q, mus);

    /* Sample */
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      auto it_inter = inter_output.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec new_state = mvrnorm(it_inter->mu, it_inter->Sigma_chol);
      ans.new_particle(std::move(new_state), &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.get_state(), it_inter->mu, it_inter->sigma_chol_inv);

      debug_msg_while_sampling(data, p, it_inter->mu, it_inter->Sigma_chol);
    }

    return(ans);
  }
};


#endif
