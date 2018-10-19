#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

#include "PF_utils.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"


/* sample function */
struct sample_output {
  arma::vec sample;
  double log_importance_dens;
};

class proposal_sampler {
public:
  virtual sample_output
  operator()(const arma::mat&, const arma::mat&) const = 0;

  virtual ~proposal_sampler() = default;
};

class proposal_sampler_mvtrnorm : public proposal_sampler {
  const int nu;

public:
  proposal_sampler_mvtrnorm() = delete;
  proposal_sampler_mvtrnorm(const int);

  sample_output operator()(const arma::mat&, const arma::mat&) const override;
};

class proposal_sampler_mvnrnorm : public proposal_sampler {
public:
  sample_output operator()(const arma::mat&, const arma::mat&) const override;
};

std::unique_ptr<proposal_sampler> get_sampler(const PF_data &data);

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
    const arma::mat *Q_chol, *Q_chol_inv;
    const arma::vec *mean;

    if(is_forward){
      Q_chol     = &data.Q_0.chol();
      Q_chol_inv = &data.Q_0.chol_inv();
      mean = &data.a_0;

    } else {
      Q_chol     = &data.uncond_covar_state(data.d + 1).chol();
      Q_chol_inv = &data.uncond_covar_state(data.d + 1).chol_inv();
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

    if(data.nu > 0L){
      proposal_sampler_mvtrnorm sampler(data.nu);

      double max_weight = -std::numeric_limits<double>::max();
      for(arma::uword i = 0; i < data.N_first; ++i){
        auto smp = sampler(*Q_chol, *Q_chol_inv);
        ans.new_particle(smp.sample + *mean, nullptr);
        ans[i].log_weight =
          dmvnrm_log(smp.sample, *Q_chol_inv) - smp.log_importance_dens;
        max_weight = std::max(max_weight, smp.log_importance_dens);
      }

      normalize_log_weights<false, true>(ans, max_weight);

    } else {
      double log_weight = log(1. / data.N_first);
      for(arma::uword i = 0; i < data.N_first; ++i){
        arma::vec err = mvrnorm(*Q_chol);
        ans.new_particle(err + *mean, nullptr);
        ans[i].log_weight = log_weight;
      }

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

    std::unique_ptr<covarmat> Q_use;
    if(is_forward)
      Q_use.reset(new covarmat(
          data.Q.mat()           + data.Q_proposal_xtra.mat()));
    else
      Q_use.reset(new covarmat(
          data.bw_covar(t).mat() + data.Q_proposal_xtra.mat()));

    if(data.debug > 2)
      data.log(3)
        << "Sampling new cloud from normal distribution with chol(Q) given by"
        << std::endl << Q_use->chol();

    auto it = resample_idx.begin();
    std::unique_ptr<proposal_sampler> sampler = get_sampler(data);
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      arma::vec mu;
      if(is_forward){
        mu = data.state_trans->map(cl[*it].get_state()).sv;

      } else {
        mu = data.bw_mean(t, cl[*it].get_state());

      }

      auto smp = (*sampler)(Q_use->chol(), Q_use->chol_inv());
      ans.new_particle(data.err_state->map(smp.sample).sv + mu, &cl[*it]);
      particle &p = ans[i];
      p.log_importance_dens = smp.log_importance_dens;

    }

    return ans;
  }

  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS){
    cloud ans;
    ans.reserve(data.N_smooth);

    bw_fw_particle_combiner combiner(data);

    covarmat rng_mat(combiner.Q.mat() + data.Q_proposal_xtra_state.mat());

    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by" << std::endl
                  << rng_mat.chol();
    }

    std::unique_ptr<proposal_sampler> sampler = get_sampler(data);
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mu = combiner(fw_p, bw_p);
      auto smp = (*sampler)(rng_mat.chol(), rng_mat.chol_inv());
      ans.new_particle(smp.sample + mu, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = smp.log_importance_dens;
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
    arma::vec parent = cl.get_weigthed_mean();

    /* compute means and covariances */
    auto inter_output = taylor_normal_approx_w_cloud_mean
      (dens_calc, data, t, data.Q, parent, cl, is_forward);

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
      rng_covar.reset(new covarmat(
          inter_output.Sigma + data.Q_proposal_xtra.mat()));
    else
      rng_covar.reset(new covarmat(
        inter_output.Sigma + data.Q_proposal_xtra_state.mat()));
    debug_msg_before_sampling(data, rng_covar->chol(), inter_output.mu);

    std::unique_ptr<proposal_sampler> sampler = get_sampler(data);
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      auto it = resample_idx.begin() + i;
      arma::vec &mu_j = inter_output.mu_js[*it];

      auto smp = (*sampler)(rng_covar->chol(), rng_covar->chol_inv());
      arma::vec &err = smp.sample;
      if(is_forward)
        ans.new_particle(mu_j + data.err_state->map(err).sv, &cl[*it]);
      else
        ans.new_particle(mu_j +                     err    , &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = smp.log_importance_dens;

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }

  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS){
    /* Find weighted mean estimate */
    bw_fw_particle_combiner combiner(data);

    const arma::vec fw_mean = fw_cloud.get_weigthed_mean(),
                    bw_mean = bw_cloud.get_weigthed_mean();

    const arma::vec alpha_bar = combiner(fw_mean, bw_mean);
    arma::vec mean_term = combiner(fw_mean, bw_mean, false);

    /* compute parts of the terms for the mean and covariance */
    auto inter_output =
      taylor_normal_approx(
          dens_calc, data, t, combiner.Q.inv(), alpha_bar, mean_term,
          2, true, false);

    /* Sample */
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    covarmat rng_covar(
        inter_output.Sigma + data.Q_proposal_xtra_state.mat());
    debug_msg_before_sampling(data, rng_covar.chol(), inter_output.mu);

    std::unique_ptr<proposal_sampler> sampler = get_sampler(data);
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      auto it_fw = fw_idx.begin() + i;
      auto it_bw = bw_idx.begin() + i;
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mu_j = combiner(fw_p.get_state(), bw_p.get_state(), false);

      /* add the term from the Taylor approximation and draw error term */
      mu_j = solve_w_precomputed_chol(inter_output.Sigma_inv_chol, mu_j) +
        inter_output.mu;

      auto smp = (*sampler)(rng_covar.chol(), rng_covar.chol_inv());
      ans.new_particle(smp.sample + mu_j, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = smp.log_importance_dens;

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
      const PF_data &data, const particle &p, const arma::vec &mu,
      const arma::mat &Sigma_chol){
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
    std::unique_ptr<proposal_sampler> sampler = get_sampler(data);
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      auto it = resample_idx.begin() + i;
      auto &inter_o = inter_output[*it];

      if(is_forward)
        rng_covar.reset(new covarmat(
            inter_o.sigma + data.Q_proposal_xtra.mat()));
      else
        rng_covar.reset(new covarmat(
            inter_o.sigma + data.Q_proposal_xtra_state.mat()));

      auto smp = (*sampler)(rng_covar->chol(), rng_covar->chol_inv());
      arma::vec &err = smp.sample;
      if(is_forward)
        ans.new_particle(data.err_state->map(err).sv + inter_o.mu, &cl[*it]);
      else
        ans.new_particle(                    err     + inter_o.mu, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = smp.log_importance_dens;

      debug_msg_while_sampling(data, p, inter_o.mu, rng_covar->chol());
    }

    return(ans);
  }

  static cloud sample_smooth(SAMPLE_SMOOTH_ARGS){
    bw_fw_particle_combiner combiner(data);

    auto begin_fw = fw_idx.begin();
    auto begin_bw = bw_idx.begin();

    cloud ans;
    ans.reserve(data.N_smooth);
    std::vector<std::unique_ptr<covarmat>> rng_mats(data.N_smooth);
    std::vector<arma::vec> propsal_means(data.N_smooth);

    /* find input for sampling */
    arma::uvec r_set = get_risk_set(data, t);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      const particle &fw_p = fw_cloud[*(begin_fw + i)];
      const particle &bw_p = bw_cloud[*(begin_bw + i)];

      const arma::vec alpha_bar = combiner(fw_p, bw_p);
      arma::vec mean_term = combiner(fw_p, bw_p, false);

      /* compute parts of the terms for the mean and covariance */
      auto inter_output =
      taylor_normal_approx(
        dens_calc, data, t, combiner.Q.inv(), alpha_bar, mean_term, r_set,
        5, false, false);

      std::unique_ptr<covarmat> new_ptr(
          new covarmat(inter_output.Sigma + data.Q_proposal_xtra_state.mat()));
      std::swap(rng_mats[i], new_ptr);
      propsal_means[i] =
        solve_w_precomputed_chol(inter_output.Sigma_inv_chol, mean_term) +
        inter_output.mu;
    }

    /* Sample */
    std::unique_ptr<proposal_sampler> sampler = get_sampler(data);
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      const particle &fw_p = fw_cloud[*(begin_fw + i)];
      const particle &bw_p = bw_cloud[*(begin_bw + i)];

      auto smp = (*sampler)(rng_mats[i]->chol(), rng_mats[i]->chol_inv());
      ans.new_particle(smp.sample + propsal_means[i], &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = smp.log_importance_dens;

      debug_msg_while_sampling(data, p, propsal_means[i], rng_mats[i]->chol());
    }

    return(ans);
  }
};

#undef SAMPLE_SMOOTH_ARGS
#undef SAMPLE_COMMON_ARGS
#endif
