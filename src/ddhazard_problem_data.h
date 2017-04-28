#ifndef DDHAZARD_PROBLEM_DATA
#define DDHAZARD_PROBLEM_DATA

#include "arma_n_rcpp.h"
#include <future>

class problem_data{
public:
  // Initalize constants
  const int d;
  const Rcpp::List risk_sets;

  const int n_params_state_vec_fixed;
  const int n_params_state_vec_varying; // NB: not including order
  const int n_params_state_vec;         // NB: not including order
  const int space_dim_in_arrays;

  const arma::span span_current_cov;         // indicies for dot product for linear predictor
  const arma::span span_current_cov_varying; // indicies for varying parameters in dot product
  const arma::span span_fixed_params;        // indicies of fixed terms in state vector

  const arma::mat &F_;
  const arma::mat T_F_;

  const bool any_dynamic;
  const bool any_fixed_in_M_step;

  arma::mat X;
  arma::mat fixed_terms; // used if fixed terms are estimated in the M-step
  const arma::vec &weights;

  const std::vector<double> I_len;
  const double event_eps; // something small

  const int n_threads;

  const arma::vec &tstart;
  const arma::vec &tstop;
  const arma::ivec &is_event_in_bin;

  const double min_start;

  const double eps_fixed_parems;
  const int max_it_fixed_params;

  const double denom_term;

  const bool debug;
  const double LR;

  const bool any_fixed_terms_in_state_vec;
  const bool use_pinv;
  const std::string criteria;

  const int EKF_batch_size;

  // Declare non constants. Some are intialize
  arma::vec fixed_parems;

  arma::mat &Q;
  arma::mat &Q_0;

  arma::mat a_t_t_s;
  arma::mat a_t_less_s;

  arma::cube V_t_t_s;
  arma::cube V_t_less_s;
  arma::cube B_s;

  arma::cube lag_one_cov;

  problem_data(const int n_fixed_terms_in_state_vec_,
               arma::mat &X_,
               arma::mat &fixed_terms_,
               const arma::vec &tstart_,
               const arma::vec &tstop_, const arma::ivec &is_event_in_bin_,
               const arma::colvec &a_0,
               const arma::vec &fixed_parems_start,
               arma::mat &Q_0_,
               arma::mat &Q_,
               const Rcpp::List &risk_obj,
               const arma::mat &F__,
               const double eps_fixed_parems_,
               const int max_it_fixed_params_,
               const arma::vec &weights_,
               const int n_max, const double eps,
               const bool verbose,
               const int order_, const bool est_Q_0,
               const bool debug_,
               Rcpp::Nullable<Rcpp::NumericVector> LR_,
               const int n_threads_,
               const double denom_term_,
               const bool use_pinv_,
               const std::string criteria_,
               const int EKF_batch_size_):
    d(Rcpp::as<int>(risk_obj["d"])),
    risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"])),

    n_params_state_vec_fixed(n_fixed_terms_in_state_vec_),
    n_params_state_vec_varying((a_0.size() - n_fixed_terms_in_state_vec_) / order_),
    n_params_state_vec(n_params_state_vec_fixed + n_params_state_vec_varying),
    space_dim_in_arrays(n_params_state_vec_varying * order_ + n_params_state_vec_fixed),

    span_current_cov(0, n_params_state_vec - 1),
    span_current_cov_varying(0, n_params_state_vec_varying - 1),
    span_fixed_params(n_params_state_vec_varying, n_params_state_vec - 1),

    F_(F__),
    T_F_(F_.t()),

    any_dynamic(X_.n_elem > 0),
    any_fixed_in_M_step(fixed_terms_.n_elem > 0),

    X(X_.begin(), X_.n_rows, X_.n_cols, false),
    fixed_terms(fixed_terms_.begin(), fixed_terms_.n_rows, fixed_terms_.n_cols, false),
    weights(weights_),
    I_len(Rcpp::as<std::vector<double> >(risk_obj["I_len"])),

    event_eps(d * std::numeric_limits<double>::epsilon()),

    n_threads((n_threads_ > 0) ? n_threads_ : std::thread::hardware_concurrency()),

    tstart(tstart_),
    tstop(tstop_),
    is_event_in_bin(is_event_in_bin_),
    min_start(Rcpp::as<double>(risk_obj["min_start"])),

    eps_fixed_parems(eps_fixed_parems_),
    max_it_fixed_params(max_it_fixed_params_),
    denom_term(denom_term_),

    debug(debug_),
    LR(LR_.isNotNull() ? Rcpp::as< Rcpp::NumericVector >(LR_)[0] : 1.0),

    any_fixed_terms_in_state_vec(n_fixed_terms_in_state_vec_ > 0),
    use_pinv(use_pinv_),
    criteria(criteria_),

    EKF_batch_size(EKF_batch_size_),

    fixed_parems(fixed_parems_start.begin(), fixed_parems_start.n_elem),
    Q(Q_),
    Q_0(Q_0_)
  {
    a_t_t_s = arma::mat(space_dim_in_arrays, d + 1);
    a_t_less_s = arma::mat(space_dim_in_arrays, d);
    V_t_t_s = arma::cube(space_dim_in_arrays, space_dim_in_arrays, d + 1);
    V_t_less_s = arma::cube(space_dim_in_arrays, space_dim_in_arrays, d);
    B_s = arma::cube(space_dim_in_arrays, space_dim_in_arrays, d);

    a_t_t_s.col(0) = a_0;

    lag_one_cov = arma::cube(space_dim_in_arrays, space_dim_in_arrays, d);

    if(debug)
      Rcpp::Rcout << "Using " << n_threads << " threads" << std::endl;
  }

  problem_data & operator=(const problem_data&) = delete;
  problem_data(const problem_data&) = delete;
  problem_data() = delete;
};

// Class with further members for the Extended Kalman Filter
class problem_data_EKF : public problem_data{
public:
  // locks for parallel implementation when we need to perform reduction from
  // local score vectors and information matricies
  std::mutex m_u;
  std::mutex m_U;

  const int n_in_last_set;
  const bool is_mult_NR;
  const double NR_eps;
  const unsigned int NR_it_max;

  // Vector for score and information matrix
  arma::colvec u;
  arma::mat U;

  // Needed for lag one covariance
  arma::mat z_dot;
  arma::vec H_diag_inv;
  arma::mat K_d;

  problem_data_EKF(const int n_fixed_terms_in_state_vec_,
                   arma::mat &X, arma::mat &fixed_terms,
                   const arma::vec &tstart_,
                   const arma::vec &tstop_, const arma::ivec &is_event_in_bin_,
                   const arma::colvec &a_0,
                   const arma::vec &fixed_parems_start,
                   arma::mat &Q_0_,
                   arma::mat &Q_,
                   const Rcpp::List &risk_obj,
                   const arma::mat &F__,
                   Rcpp::Nullable<Rcpp::NumericVector> NR_eps_,
                   Rcpp::Nullable<Rcpp::NumericVector> LR_,
                   const double eps_fixed_parems_,
                   const int max_it_fixed_params_,
                   const arma::vec &weights_,
                   const int n_max, const double eps,
                   const bool verbose,
                   const int order_, const bool est_Q_0,
                   const bool is_cont_time,
                   const unsigned int NR_it_max_,
                   const bool debug_,
                   const int n_threads_,
                   const double denom_term_,
                   const bool use_pinv_,
                   const std::string criteria_,
                   const int EKF_batch_size_):
    problem_data(n_fixed_terms_in_state_vec_, X, fixed_terms, tstart_, tstop_, is_event_in_bin_, a_0,
                 fixed_parems_start, Q_0_, Q_, risk_obj, F__,
                 eps_fixed_parems_, max_it_fixed_params_,
                 weights_,
                 n_max, eps, verbose,
                 order_, est_Q_0, debug_,
                 LR_,
                 n_threads_, denom_term_,
                 use_pinv_, criteria_, EKF_batch_size_),

                 n_in_last_set(Rcpp::as<arma::uvec>(risk_sets[d - 1]).size()),
                 is_mult_NR(NR_eps_.isNotNull()),
                 NR_eps(is_mult_NR ? Rcpp::as< Rcpp::NumericVector >(NR_eps_)[0] : 0.0),
                 NR_it_max(NR_it_max_)
  {
    u = arma::colvec(space_dim_in_arrays);
    U = arma::mat(space_dim_in_arrays, space_dim_in_arrays);

    z_dot = arma::mat(space_dim_in_arrays, n_in_last_set * (is_cont_time + 1));
    H_diag_inv = arma::vec(n_in_last_set * (is_cont_time + 1));
  }

  problem_data_EKF & operator=(const problem_data_EKF&) = delete;
  problem_data_EKF(const problem_data_EKF&) = delete;
  problem_data_EKF() = delete;
};

#endif
