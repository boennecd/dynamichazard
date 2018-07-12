#include "PF_utils.h"


inline covarmat set_bw_fw_particle_combiner_Q
  (const PF_data &data, const covarmat &Q_trans){
  arma::mat Q_inv_arg = data.err_state_inv->map(
    Q_trans.inv, both, trans).sv;
  Q_inv_arg += data.state_trans->map(Q_inv_arg, both, trans).sv;
  return covarmat(Q_inv_arg.i());
}

bw_fw_particle_combiner::bw_fw_particle_combiner
  (const PF_data &data):
  data(data), Q_trans(data.Q_proposal_smooth),
  Q(set_bw_fw_particle_combiner_Q(data, Q_trans)) {}

arma::vec bw_fw_particle_combiner::operator()
  (const particle &fw_p, const particle &bw_p) const {
  return this->operator()(fw_p.get_state(), bw_p.get_state());
}

arma::vec bw_fw_particle_combiner::operator()
  (const arma::vec &fw_p, const arma::vec &bw_p) const {
  /* compute part of the mean from the forward particle */
  arma::vec mu = data.state_trans_err->map(fw_p).sv;
  mu = solve_w_precomputed_chol(Q_trans.chol, mu);
  mu = data.state_trans_err->map(mu, trans).sv;

  /* add part of the mean from the backward particle */
  {
    arma::vec bw_term = data.err_state_inv->map(bw_p).sv;
    bw_term = solve_w_precomputed_chol(Q_trans.chol, bw_term);
    mu += data.state_trans_err->map(bw_term, trans).sv;
  }

  return Q.mat * mu;
}
