#include "importance_samplers.h"

std::unique_ptr<proposal_sampler> get_sampler(const PF_data &data){
  if(data.nu > 0L)
    return std::unique_ptr<proposal_sampler>(
      new proposal_sampler_mvtrnorm(data.nu));

  return
    std::unique_ptr<proposal_sampler>(new proposal_sampler_mvnrnorm());
}




proposal_sampler_mvtrnorm::proposal_sampler_mvtrnorm(const int nu):
  nu(nu) {}

sample_output proposal_sampler_mvtrnorm::operator()(
    const arma::mat &chol, const arma::mat &chol_inv) const
{
  sample_output out;
  out.sample = mvtrnorm(chol, nu);
  out.log_importance_dens = dmvtrm_log(out.sample, chol_inv, nu);
  return out;
}



sample_output proposal_sampler_mvnrnorm::operator()(
    const arma::mat &chol, const arma::mat &chol_inv) const
{
  sample_output out;
  out.sample = mvrnorm(chol);
  out.log_importance_dens = dmvnrm_log(out.sample, chol_inv);
  return out;
}
