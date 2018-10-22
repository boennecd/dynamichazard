#include "PF_utils.h"
#include "../sample_funcs.h"

inline covarmat set_bw_fw_particle_combiner_Q
  (const PF_data &data, const covarmat &Q_trans){
  arma::mat Q_inv_arg = data.err_state->map(Q_trans.inv()).sv;
  Q_inv_arg += data.state_trans->map(Q_inv_arg, both, trans).sv;
  return covarmat(Q_inv_arg.i());
}

bw_fw_particle_combiner::bw_fw_particle_combiner
  (const PF_data &data):
  data(data), Q_trans(data.Q),
  Q(set_bw_fw_particle_combiner_Q(data, Q_trans)) {}

arma::vec bw_fw_particle_combiner::operator()
  (const particle &fw_p, const particle &bw_p, const bool do_transform) const {
  return this->operator()(fw_p.get_state(), bw_p.get_state(), do_transform);
}

arma::vec bw_fw_particle_combiner::operator()
  (const arma::vec &fw_p, const arma::vec &bw_p, const bool do_transform) const {
  /* compute part of the mean from the forward particle
   * R Q^{-1}R^\top Fx_{t - 1}
   * TODO: issues w/ higher orders currently... */
  arma::vec mu = data.state_trans_err->map(fw_p).sv;
  mu = solve_w_precomputed_chol(Q_trans.chol(), mu);
  mu = data.err_state->map(mu).sv;

  /* add part of the mean from the backward particle */
  {
    /* F^\top R Q^{-1} R^\top x_{t + 1}
     * TODO: issues w/ higher orders currently... */
    arma::vec bw_term = data.err_state_inv->map(bw_p).sv;
    bw_term = solve_w_precomputed_chol(Q_trans.chol(), bw_term);
    mu += data.state_trans_err->map(bw_term, trans).sv;
  }

  if(!do_transform)
    return mu;

  return Q.mat() * mu;
}


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
