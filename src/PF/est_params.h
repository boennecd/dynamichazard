#ifndef PF_EST_PARAMS
#define PF_EST_PARAMS
#include "PF_utils.h"

struct PF_summary_stats {
  /* E[X_t] */
  std::vector<arma::vec> E_xs;
  /* E[(X_t - FX_{t-1})(X_t - FX_{t-1})^T] */
  std::vector<arma::mat> E_x_less_x_less_one_outers;
};

PF_summary_stats
  compute_PF_summary_stats
  (const smoother_output&, const arma::vec&, const arma::mat&,
   const arma::mat&, const arma::mat F = arma::mat(),
   const bool do_use_F = false, const bool do_compute_E_x = true);

struct PF_parameters {
  /* time zero value of the space vector */
  arma::mat a_0;
  /* R^\top F matrix */
  arma::mat R_top_F;
  /* covariance matrix in state equation */
  arma::mat Q;
  /* output form QR. Names are confusing though */
  arma::mat R;
  arma::mat F;
  arma::mat dev;
};

PF_parameters
  est_params_dens
  (const smoother_output&, const arma::vec&, const arma::mat&,
   const arma::mat&, const arma::mat&, const int, const bool, const bool,
   const unsigned long max_bytes = 5000000,
   const bool only_QR = false);

#endif
