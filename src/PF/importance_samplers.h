#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

#include "PF_utils.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"


#define SAMPLE_SMOOTH_ARGS                                \
  pf_base_dens &dens_calc,                                \
  const PF_data &data,                                    \
  cloud &fw_cloud,                                        \
  const arma::uvec &fw_idx,                               \
  cloud &bw_cloud,                                        \
  const arma::uvec &bw_idx,                               \
  const unsigned int t

#define SAMPLE_COMMON_ARGS                                     \
  pf_base_dens &dens_calc,                                     \
  const PF_data &data,                                         \
  cloud &cl,                                                   \
  const arma::uvec &resample_idx,                              \
  const unsigned int t

/*
  Each "importance_sampler" has two static functions:
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
  static cloud sample_first_state_n_set_weights
  (pf_base_dens &dens_calc, const PF_data &data){
    cloud ans;
    ans.reserve(data.N_first);
    const arma::mat *Q_chol;
    const arma::vec *mean;

    if(is_forward){
      Q_chol = &data.Q_0.chol();
      mean = &data.a_0;

    } else {
      Q_chol = &data.uncond_covar_state(data.d + 1).chol();
      mean = &data.uncond_mean_state(data.d + 1);

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
      arma::vec err = mvrnorm(*Q_chol);
      ans.new_particle(err + *mean, nullptr);
      ans[i].log_weight = log_weight;
    }

    return(ans);
  }
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
  static cloud sample(SAMPLE_COMMON_ARGS, nothing unused){
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    const covarmat *Q_use;
    if(is_forward){
      Q_use = &data.Q_proposal;

    } else {
      Q_use = &data.bw_covar(t);

    }

    if(data.debug > 2){
      data.log(3)
        << "Sampling new cloud from normal distribution with chol(Q) given by"
        << std::endl << Q_use->chol();

    }

    auto it = resample_idx.begin();
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      arma::vec mu;
      if(is_forward){
        mu = data.state_trans->map(cl[*it].get_state()).sv;

      } else {
        mu = data.bw_mean(t, cl[*it].get_state());

      }

      arma::vec err = mvrnorm(Q_use->chol());
      ans.new_particle(data.err_state->map(err).sv + mu, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(err, Q_use->chol_inv());

    }

    return ans;
  }

  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS){
    cloud ans;
    ans.reserve(data.N_smooth);

    bw_fw_particle_combiner combiner(data);

    covarmat rng_mat(combiner.Q.mat() + data.Q_proposal_state.mat());

    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by" << std::endl
                  << rng_mat.chol();
    }

    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mu = combiner(fw_p, bw_p);
      arma::vec err = mvrnorm(rng_mat.chol());
      ans.new_particle(err + mu, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(err, rng_mat.chol_inv());
    }

    return ans;
  }
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
  inline static void debug_msg_before_sampling(
      const PF_data &data,
      const arma::mat &Sigma_chol, const arma::vec mu){
    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Sigma) given by" << std::endl
                  << Sigma_chol
                  << "The mean before accounting for the parent (and child) particle is:" << std::endl
                  << mu.t();
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
  static cloud sample(SAMPLE_COMMON_ARGS, nothing unused){
    /* Find weighted mean estimate */
    arma::vec alpha_bar = cl.get_weigthed_mean();

    if(is_forward)
      alpha_bar = data.state_trans->map(alpha_bar).sv;
    else
      alpha_bar = data.bw_mean(t, alpha_bar);

    /* compute means and covariances */
    auto inter_output = taylor_normal_approx_w_cloud_mean
      (dens_calc, data, t, data.Q, alpha_bar, cl, is_forward);

    return(sample(dens_calc, data, cl, resample_idx, t, inter_output));
  }

  static cloud sample(
      SAMPLE_COMMON_ARGS,
      input_for_normal_apprx_w_cloud_mean &inter_output){
    /* Sample */
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    std::unique_ptr<covarmat> rng_covar;
    if(is_forward)
      rng_covar.reset(new covarmat(inter_output.Sigma + data.Q_proposal.mat()));
    else
      rng_covar.reset(new covarmat(
        inter_output.Sigma + data.Q_proposal_state.mat()));
    debug_msg_before_sampling(data, rng_covar->chol(), inter_output.mu);

    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      auto it = resample_idx.begin() + i;
      arma::vec &mu_j = inter_output.mu_js[*it];

      arma::vec err = mvrnorm(rng_covar->chol());
      if(is_forward)
        ans.new_particle(mu_j + data.err_state->map(err).sv, &cl[*it]);
      else
        ans.new_particle(mu_j +                     err    , &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(err, rng_covar->chol_inv());

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }

  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS){
    /* Find weighted mean estimate */
    bw_fw_particle_combiner combiner(data);

    const arma::vec alpha_bar = combiner(
      fw_cloud.get_weigthed_mean(), bw_cloud.get_weigthed_mean());

    /* compute parts of the terms for the mean and covariance */
    auto inter_output =
      taylor_normal_approx(
          dens_calc, data, t, combiner.Q.inv(), alpha_bar, 2, true, false);

    /* Sample */

    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    covarmat rng_covar(
        inter_output.Sigma + data.Q_proposal_state.mat());
    debug_msg_before_sampling(data, rng_covar.chol(), inter_output.mu);
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      /* compute part of the mean from the forward particle */
      arma::vec mu_j = data.state_trans_err->map(fw_p.get_state()).sv;
      mu_j = solve_w_precomputed_chol(combiner.Q_trans.chol(), mu_j);
      mu_j = data.state_trans_err->map(mu_j, trans).sv;

      /* add part of the mean from the backward particle */
      {
        arma::vec bw_term = data.err_state_inv  ->map(bw_p.get_state()).sv;
        bw_term = solve_w_precomputed_chol(combiner.Q_trans.chol(), bw_term);
        mu_j += data.state_trans_err            ->map(bw_term, trans).sv;
      }

      /* add the term from the Taylor approximation and draw error term */
      mu_j = solve_w_precomputed_chol(inter_output.Sigma_inv_chol, mu_j) +
        inter_output.mu;

      arma::vec err = mvrnorm(rng_covar.chol());
      ans.new_particle(err + mu_j, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens =
        dmvnrm_log(err, rng_covar.chol_inv());

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }
};

/*
  Sampler which makes a normal approximation for the observed outcome made
  around the mean of the parent particle.
*/

template<bool is_forward>
class importance_dens_normal_approx_w_particles  :
  public importance_dens_base<is_forward> {

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
  static cloud sample(SAMPLE_COMMON_ARGS, nothing unused) {
    /* compute means and covariances */
    auto inter_output =
      taylor_normal_approx_w_particles
      (dens_calc, data, t, data.Q, cl, is_forward);

    return(sample(dens_calc, data, cl, resample_idx, t, inter_output));
  }

  static cloud sample(
      SAMPLE_COMMON_ARGS,
      input_for_normal_apprx_w_particle_mean &inter_output){
    /* Sample */
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    std::unique_ptr<covarmat> rng_covar;
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      auto it = resample_idx.begin() + i;
      auto &inter_o = inter_output[*it];

      if(is_forward)
        rng_covar.reset(new covarmat(inter_o.sigma + data.Q_proposal.mat()));
      else
        rng_covar.reset(new covarmat(
            inter_o.sigma + data.Q_proposal_state.mat()));

      arma::vec err = mvrnorm(rng_covar->chol());
      if(is_forward)
        ans.new_particle(data.err_state->map(err).sv + inter_o.mu, &cl[*it]);
      else
        ans.new_particle(                    err     + inter_o.mu, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(err, rng_covar->chol_inv());

      debug_msg_while_sampling(data, p, inter_o.mu, rng_covar->chol());
    }

    return(ans);
  }

  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS){
    /* Compute means before accounting for outcomes */
    std::vector<arma::vec> mus(data.N_smooth);

    bw_fw_particle_combiner combiner(data);

    auto begin_fw = fw_idx.begin();
    auto begin_bw = bw_idx.begin();
    auto begin_mu = mus.begin();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      arma::vec &mu_j = *(begin_mu + i);
      mu_j = combiner(fw_cloud[*(begin_fw + i)], bw_cloud[*(begin_bw + i)]);
    }

    /* compute means and covariances */
    auto inter_output =
      taylor_normal_approx_w_particles
      (dens_calc, data, t, combiner.Q, mus, false);

    /* Sample */
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      auto it_inter = inter_output.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      covarmat rng_mat(it_inter->sigma + data.Q_proposal_state.mat());

      arma::vec err = mvrnorm(rng_mat.chol());
      ans.new_particle(err + it_inter->mu, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(err, rng_mat.chol_inv());

      debug_msg_while_sampling(data, p, it_inter->mu, rng_mat.chol());
    }

    return(ans);
  }
};

#undef SAMPLE_SMOOTH_ARGS
#undef SAMPLE_COMMON_ARGS
#endif
