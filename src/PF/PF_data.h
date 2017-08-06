#ifndef PF_DATA
#define PF_DATA

#include "../arma_n_rcpp.h"
#include "../problem_data.h"

// data holder for particle filtering
class PF_data : public problem_data {
public:
  /* Number of paprticles in forward and/or backward filter */
  const arma::uword N_fw_n_bw;
  const arma::uword N_smooth;
  const double forward_backward_ESS_threshold;

  const arma::vec &a_0;

  const arma::mat Q_chol;
  const arma::mat Q_half_chol;
  const arma::mat Q_chol_inv;
  const arma::mat Q_half_chol_inv;

  PF_data(const int n_fixed_terms_in_state_vec_,
          arma::mat &X_,
          arma::mat &fixed_terms_,
          const arma::vec &tstart_,
          const arma::vec &tstop_, const arma::ivec &is_event_in_bin_,
          const arma::colvec &a_0,
          arma::mat &Q_0_,
          arma::mat &Q_,
          const Rcpp::List &risk_obj,
          const arma::mat &F__,
          const int n_max,
          const int order_,
          const int n_threads_,

          // new arguments
          const arma::uword N_fw_n_bw,
          const arma::uword N_smooth,
          const double forward_backward_ESS_threshold):
    problem_data(
      n_fixed_terms_in_state_vec_,
      X_,
      fixed_terms_,
      tstart_,
      tstop_, is_event_in_bin_,
      a_0,
      Q_0_,
      Q_,
      risk_obj,
      F__,
      n_max,
      order_,
      n_threads_),

      N_fw_n_bw(N_fw_n_bw),
      N_smooth(N_smooth),
      forward_backward_ESS_threshold(forward_backward_ESS_threshold),

      a_0(a_0),

      Q_chol(arma::chol(Q)),
      Q_half_chol(Q_half_chol / sqrt(2.)),
      Q_chol_inv(arma::inv(arma::trimatu(Q_chol))),
      Q_half_chol_inv(arma::inv(arma::trimatu(Q_half_chol)))
    {}

  PF_data & operator=(const PF_data&) = delete;
  PF_data(const PF_data&) = delete;
  PF_data() = delete;
};

#endif
