#include "PF_utils.h"
#include "../sample_funcs.h"

cloud re_sample_cloud(const unsigned int size, const cloud cl){
  if(size >= cl.size())
    Rcpp::stop("size greater than or equal to cl.size() in 're_sample_cloud'");

  arma::vec probs(cl.size());
  double *p = probs.begin();
  for(auto it = cl.begin(); it != cl.end(); ++it, ++p)
    *p = std::exp(it->log_weight);

  std::map<arma::uword, arma::uword> idx =
    sample_n_count_replicas<systematic_resampling>(size, probs);

  cloud out;
  out.reserve(idx.size());
  unsigned int i = 0;
  for (auto it = idx.begin(); it != idx.end(); it++, i++)
  {
    const particle &to_copy = cl[it->first];
    out.new_particle(to_copy.get_state(), to_copy.parent, to_copy.child);
    particle &p = out[i];
    p.log_importance_dens = to_copy.log_importance_dens;
    p.log_likelihood_term = to_copy.log_likelihood_term;
    p.log_weight = log(((double)it->second) / size);
  }

  return out;
}



template<bool is_forward>
get_approx_use_mean_output get_approx_use_mean(
    std::shared_ptr<PF_cdist> y_dist, cloud &PF_cloud, const PF_data &data,
    pf_dens &dens_calc, arma::uword t){
  unsigned int n_elem = PF_cloud.size();
  std::vector<std::unique_ptr<dist_comb>> out(n_elem);

  std::unique_ptr<PF_cdist> other;
  std::shared_ptr<PF_cdist> prior;
  arma::vec other_state = PF_cloud.get_weigthed_mean(), start;
  std::vector<PF_cdist*> objs;
  if(is_forward){
    other = dens_calc.get_fw_dist(other_state);
    start = other->get_mean();
    objs = { y_dist.get(), other.get() };

  }
  else
  {
    other = dens_calc.get_bw_dist(other_state);
    prior = dens_calc.get_prior(t);
    std::vector<PF_cdist*> start_objs = { other.get(), prior.get() };
    start = cdist_comb_generator(start_objs).get_dist_comb
      ({ &other_state })->get_mean();
    objs = { y_dist.get(), other.get(), prior.get() };

  }

  cdist_comb_generator combi_gen(
      objs, start, data.nu, &data.xtra_covar, data.covar_fac, data.ftol_rel);
  for(unsigned int i = 0; i < n_elem; ++i){ // loop over cloud elements
    auto it_cl = PF_cloud.begin() + i;
    auto it_dc = out.begin() + i;

    *it_dc = combi_gen.get_dist_comb({ (arma::vec*)&it_cl->get_state() });
  }

  return { std::move(out) };
}

template get_approx_use_mean_output get_approx_use_mean<true>(
    std::shared_ptr<PF_cdist>, cloud&, const PF_data&, pf_dens&, arma::uword);
template get_approx_use_mean_output get_approx_use_mean<false>(
    std::shared_ptr<PF_cdist>, cloud&, const PF_data&, pf_dens&, arma::uword);

template<bool is_forward>
get_approx_use_particle_output get_approx_use_particle(
    std::shared_ptr<PF_cdist> y_dist, cloud &PF_cloud, const PF_data &data,
    pf_dens &dens_calc, arma::uword t){
  unsigned int n_elem = PF_cloud.size();

  get_approx_use_particle_output out;
  out.dists = std::vector<std::unique_ptr<dist_comb>>(n_elem);

  std::unique_ptr<PF_cdist> first_dist;
  std::shared_ptr<PF_cdist> prior;
  std::vector<PF_cdist*> start_objs;
  std::unique_ptr<cdist_comb_generator> comb_start;
  if(!is_forward){
    first_dist = dens_calc.get_bw_dist(PF_cloud.begin()->get_state());
    prior = dens_calc.get_prior(t);
    start_objs = { first_dist.get(), prior.get() };
    comb_start = std::unique_ptr<cdist_comb_generator>(
      new cdist_comb_generator(start_objs));

  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned int i = 0; i < n_elem; ++i){ // loop over cloud elements
    auto it_cl = PF_cloud.begin() + i;
    auto it_dc = out.dists.begin() + i;

    std::unique_ptr<PF_cdist> other;
    std::vector<PF_cdist*> objs;
    arma::vec start;
    if(is_forward){
      other = dens_calc.get_fw_dist(it_cl->get_state());
      start = other->get_mean();
      objs = { y_dist.get(), other.get() };

    }
    else {
      other = dens_calc.get_bw_dist(it_cl->get_state());
      start = comb_start->get_dist_comb
        ({ (arma::vec*)&it_cl->get_state() })->get_mean();
      objs = { y_dist.get(), other.get(), prior.get() };

    }

    cdist_comb_generator combi_gen(
        objs, start, data.nu, &data.xtra_covar, data.covar_fac, data.ftol_rel);

    *it_dc =
      combi_gen.get_dist_comb({ (arma::vec*)&it_cl->get_state() });
  }

  return out;
}
template get_approx_use_particle_output
  get_approx_use_particle<true>(
    std::shared_ptr<PF_cdist>, cloud&, const PF_data&, pf_dens&, arma::uword);
template get_approx_use_particle_output
  get_approx_use_particle<false>(
    std::shared_ptr<PF_cdist>, cloud&, const PF_data&, pf_dens&, arma::uword);


double get_weight_from_particle(const particle &p) {
  return p.log_weight;
}
double get_resample_weight_from_particle(const particle &p){
  return p.log_resampling_weight;
}
