#ifndef PF_EST_PARAMS
#define PF_EST_PARAMS
#include "PF_utils.h"

struct PF_summary_stats_RW {
  /* E[X_t] */
  std::vector<arma::vec> E_xs;
  /* E[(X_t - FX_{t-1})(X_t - FX_{t-1})^T] */
  std::vector<arma::mat> E_x_less_x_less_one_outers;
};

PF_summary_stats_RW
  compute_summary_stats_first_o_RW
  (const smoother_output&, const arma::vec&, const arma::mat&,
   const arma::mat&);

struct PF_parameters {
  /* time zero value of the space vector */
  arma::mat a_0;
  /* R^\top F matrix */
  arma::mat R_top_F;
  /* covariance matrix in state equation */
  arma::mat Q;
};

PF_parameters
  est_params_dens
  (const smoother_output&, const arma::vec&, const arma::mat&,
   const arma::mat&, const arma::mat&, int,
   const unsigned long int max_bytes = 5000000);

#endif
