#include "importance_samplers.h"

#define SAMPLE_SMOOTH_ARGS                                  \
  std::shared_ptr<PF_cdist> y_dist, pf_dens &dens_calc,     \
  const PF_data &data, cloud &fw_cloud,                     \
  const arma::uvec &fw_idx, cloud &bw_cloud,                \
  const arma::uvec &bw_idx, const arma::uword t

#define SAMPLE_COMMON_ARGS                                     \
  std::shared_ptr<PF_cdist> y_dist, pf_dens &dens_calc,        \
  const PF_data &data, cloud &cl,                              \
  const arma::uvec &resample_idx, const arma::uword t



template<bool is_forward>
cloud importance_dens_base<is_forward>::sample_first_state_n_set_weights
  (pf_dens &dens_calc, const PF_data &data){
  cloud ans;
  ans.reserve(data.N_first);

  std::unique_ptr<dist_comb> sampler, normal_dist;
  if(is_forward){
    std::shared_ptr<PF_cdist> prior = dens_calc.get_prior(0L);
    std::vector<PF_cdist*> dists = { prior.get() };
    sampler = cdist_comb_generator(dists, data.nu, &data.xtra_covar)
      .get_dist_comb({ });
    normal_dist = cdist_comb_generator(dists, -1L).get_dist_comb({ });

  } else {
    std::shared_ptr<PF_cdist> prior = dens_calc.get_prior(data.d + 1L);
    std::vector<PF_cdist*> dists = { prior.get() };
    sampler = cdist_comb_generator(dists, data.nu, &data.xtra_covar)
      .get_dist_comb({ });
    normal_dist = cdist_comb_generator(dists, -1L).get_dist_comb({ });

  }

  if(data.debug > 1){
    data.log(2) << "Sampling "
                << (is_forward ? "first" : "state d + 1")
                << " with chol(covariance matrix):" << std::endl
                << sampler->get_covar()
                << "and mean:" << std::endl
                << sampler->get_mean();
  }

  if(data.nu > 0L){
    double max_weight = -std::numeric_limits<double>::max();
    for(arma::uword i = 0; i < data.N_first; ++i){
      ans.new_particle(sampler->sample(), nullptr);
      ans[i].log_weight = normal_dist->log_density(ans[i].get_state()) -
        sampler->log_density(ans[i].get_state());

      max_weight = std::max(max_weight, ans[i].log_weight);
    }

    normalize_log_weights<false, true>(ans, max_weight);

    return ans;
  }

  double log_weight = -log(data.N_first);
  for(arma::uword i = 0; i < data.N_first; ++i){
    ans.new_particle(sampler->sample(), nullptr);
    ans[i].log_weight = log_weight;
  }

  return ans;
}

template class importance_dens_base<true>;
template class importance_dens_base<false>;



template<bool is_forward>
cloud importance_dens_no_y_dependence<is_forward>::
  sample(SAMPLE_COMMON_ARGS, nothing unused)
  {
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    std::shared_ptr<PF_cdist> d1, prior;
    std::vector<PF_cdist*> dists;
    if(is_forward){
      d1 = dens_calc.get_fw_dist(
        /* does not matter which one we use */
        cl.begin()->get_state());
      dists = { d1.get() };

    } else {
      d1 = dens_calc.get_bw_dist(
        /* does not matter which one we use */
        cl.begin()->get_state());
      prior = dens_calc.get_prior(t);
      dists = { d1.get(), prior.get() };

    }
    cdist_comb_generator comb_gen(dists, data.nu, &data.xtra_covar);

    auto it = resample_idx.begin();
    std::unique_ptr<dist_comb> sampler;
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      sampler = comb_gen.get_dist_comb({ (arma::vec*)&cl[*it].get_state() });
      ans.new_particle(sampler->sample(), &cl[*it]);

      ans[i].log_importance_dens = sampler->log_density(ans[i].get_state());
    }

    return ans;
  }

template<bool is_forward>
cloud importance_dens_no_y_dependence<is_forward>::
  sample_smooth(SAMPLE_SMOOTH_ARGS)
  {
    cloud ans;
    ans.reserve(data.N_smooth);

    std::shared_ptr<PF_cdist>
      fw_dist = dens_calc.get_fw_dist(
        /* does not matter which one we use */
        fw_cloud.begin()->get_state()),
      bw_dist = dens_calc.get_bw_dist(
        /* does not matter which one we use */
        bw_cloud.begin()->get_state());
    std::vector<PF_cdist*> dists = { fw_dist.get(), bw_dist.get() };
    cdist_comb_generator comb_gen(dists, data.nu, &data.xtra_covar);

    std::unique_ptr<dist_comb> sampler;
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      const particle &fw_p = fw_cloud[*(fw_idx.begin() + i)];
      const particle &bw_p = bw_cloud[*(bw_idx.begin() + i)];

      sampler = comb_gen.get_dist_comb
        ({ (arma::vec*)&fw_p.get_state(), (arma::vec*)&bw_p.get_state()});
      ans.new_particle(sampler->sample(), &fw_p, &bw_p);

      ans[i].log_importance_dens = sampler->log_density(ans[i].get_state());
    }

    return ans;
  }

template class importance_dens_no_y_dependence<true>;
template class importance_dens_no_y_dependence<false>;


template<bool is_forward>
void importance_dens_normal_approx_w_cloud_mean<is_forward>::
  debug_msg_while_sampling
  (const PF_data &data, const particle &p, const arma::vec &mu){
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

template<bool is_forward>
cloud importance_dens_normal_approx_w_cloud_mean<is_forward>::sample
  (SAMPLE_COMMON_ARGS, nothing unused)
  {
    std::vector<std::unique_ptr<dist_comb>> samplers =
      get_approx_use_mean<is_forward>(y_dist, cl, data, dens_calc, t);

    return sample(y_dist, dens_calc, data, cl, resample_idx, t, samplers);
  }

template<bool is_forward>
cloud importance_dens_normal_approx_w_cloud_mean<is_forward>::sample
  (SAMPLE_COMMON_ARGS, std::vector<std::unique_ptr<dist_comb>> &samplers)
  {
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
      auto it = resample_idx.begin() + i;
      std::unique_ptr<dist_comb> &sampler = samplers[*it];

      ans.new_particle(sampler->sample(), &cl[*it]);
      ans[i].log_importance_dens = sampler->log_density(ans[i].get_state());
    }

    return ans;
  }

template<bool is_forward>
cloud importance_dens_normal_approx_w_cloud_mean<is_forward>::
  sample_smooth(SAMPLE_SMOOTH_ARGS){
    /* get object to form proposal distribution */
    arma::vec fw_mean = fw_cloud.get_weigthed_mean(),
      bw_mean = bw_cloud.get_weigthed_mean();

    std::unique_ptr<PF_cdist>
      fw_dist = dens_calc.get_fw_dist(fw_mean),
      bw_dist = dens_calc.get_bw_dist(bw_mean);
    std::vector<PF_cdist*> objs =
      { y_dist.get(), fw_dist.get(), bw_dist.get() };
    arma::vec start;
    {
      std::vector<PF_cdist*> start_objs = { fw_dist.get(), bw_dist.get() };
      start = cdist_comb_generator(start_objs)
        .get_dist_comb({ &fw_mean, &bw_mean })->get_mean();
    }
    cdist_comb_generator combi_gen(objs, start, data.nu, &data.xtra_covar);

    /* get proposal distributions */
    auto begin_fw = fw_idx.begin();
    auto begin_bw = bw_idx.begin();
    std::vector<std::unique_ptr<dist_comb>> samplers(data.N_smooth);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      const particle &fw_p = fw_cloud[*(begin_fw + i)];
      const particle &bw_p = bw_cloud[*(begin_bw + i)];

      samplers[i] =
        combi_gen.get_dist_comb
        ({ (arma::vec*)&fw_p.get_state(), (arma::vec*)&bw_p.get_state() });
    }

    /* sample */
    cloud ans;
    ans.reserve(data.N_smooth);
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      const particle &fw_p = fw_cloud[*(begin_fw + i)];
      const particle &bw_p = bw_cloud[*(begin_bw + i)];
      std::unique_ptr<dist_comb> &sampler = samplers[i];

      ans.new_particle(sampler->sample(), &fw_p, &bw_p);
      ans[i].log_importance_dens = sampler->log_density(ans[i].get_state());

      if(data.debug > 4)
        debug_msg_while_sampling(data, ans[i], sampler->get_mean());
    }

    return ans;
  }

template class importance_dens_normal_approx_w_cloud_mean<true>;
template class importance_dens_normal_approx_w_cloud_mean<false>;



template<bool is_forward>
void importance_dens_normal_approx_w_particles<is_forward>::
  debug_msg_while_sampling
  (const PF_data &data, const particle &p, const arma::vec &mu,
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

template<bool is_forward>
cloud importance_dens_normal_approx_w_particles<is_forward>::sample
  (SAMPLE_COMMON_ARGS, nothing unused)
{
  std::vector<std::unique_ptr<dist_comb>> samplers =
    get_approx_use_particle<is_forward>(y_dist, cl, data, dens_calc, t);

  return sample(y_dist, dens_calc, data, cl, resample_idx, t, samplers);
}

template<bool is_forward>
cloud importance_dens_normal_approx_w_particles<is_forward>::sample
  (SAMPLE_COMMON_ARGS, std::vector<std::unique_ptr<dist_comb>> &samplers)
{
  cloud ans;
  ans.reserve(data.N_fw_n_bw);

  for(arma::uword i = 0; i < data.N_fw_n_bw; ++i){
    auto *it = resample_idx.begin() + i;
    std::unique_ptr<dist_comb> &sampler = samplers[*it];

    ans.new_particle(sampler->sample(), &cl[*it]);
    ans[i].log_importance_dens = sampler->log_density(ans[i].get_state());
  }

  return ans;
}

template<bool is_forward>
cloud importance_dens_normal_approx_w_particles<is_forward>::
  sample_smooth(SAMPLE_SMOOTH_ARGS){
    /* get object to get conditional on parent and child. It is used
     * for the starting value */
    std::unique_ptr<PF_cdist>
      fw_dum = dens_calc.get_fw_dist(fw_cloud.begin()->get_state()),
        bw_dum = dens_calc.get_bw_dist(bw_cloud.begin()->get_state());
    std::vector<PF_cdist*> start_objs = { fw_dum.get(), bw_dum.get() };
    cdist_comb_generator fw_bw_comb(start_objs);

    /* get proposal distributions */
    auto begin_fw = fw_idx.begin();
    auto begin_bw = bw_idx.begin();
    std::vector<std::unique_ptr<dist_comb>> samplers(data.N_smooth);

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      const particle &fw_p = fw_cloud[*(begin_fw + i)];
      const particle &bw_p = bw_cloud[*(begin_bw + i)];

      std::unique_ptr<PF_cdist>
        fw_dist = dens_calc.get_fw_dist(fw_p.get_state()),
          bw_dist = dens_calc.get_bw_dist(bw_p.get_state());
      std::vector<PF_cdist*> objs =
        { y_dist.get(), fw_dist.get(), bw_dist.get() };
      arma::vec start = fw_bw_comb.get_dist_comb
        ({ (arma::vec*)&fw_p.get_state(), (arma::vec*)&bw_p.get_state() })
        ->get_mean();
      cdist_comb_generator combi_gen(objs, start, data.nu, &data.xtra_covar);

      samplers[i] = combi_gen.get_dist_comb
        ({ (arma::vec*)&fw_p.get_state(), (arma::vec*)&bw_p.get_state() });
    }

    /* sample */
    cloud ans;
    ans.reserve(data.N_smooth);
    for(arma::uword i = 0; i < data.N_smooth; ++i){
      const particle &fw_p = fw_cloud[*(begin_fw + i)];
      const particle &bw_p = bw_cloud[*(begin_bw + i)];
      std::unique_ptr<dist_comb> &sampler = samplers[i];

      ans.new_particle(sampler->sample(), &fw_p, &bw_p);
      ans[i].log_importance_dens = sampler->log_density(ans[i].get_state());

      if(data.debug > 4)
        debug_msg_while_sampling(data, ans[i], sampler->get_mean(),
                                 sampler->get_covar());
    }

    return ans;
  }

template class importance_dens_normal_approx_w_particles<true>;
template class importance_dens_normal_approx_w_particles<false>;
