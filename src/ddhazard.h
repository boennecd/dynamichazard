#ifndef DDHAZARD_H
#define DDHAZARD_H

// [[Rcpp::plugins(cpp11)]]
#include "ddhazard_problem_data.h"
#include "arma_n_rcpp.h"

inline bool is_exponential_model(std::string model){
  return(model == "exp_bin" ||
         model == "exp_clip_time" ||
         model == "exp_clip_time_w_jump");
}

// Abstact solver class
class Solver {
public:
  virtual void solve() = 0;
};

// Classes for EKF method
template<typename T>
class EKF_solver : public Solver{
  problem_data_EKF &p_dat;
  const std::string model;
  unsigned long const max_threads;

  void parallel_filter_step(
      arma::uvec::const_iterator first, arma::uvec::const_iterator last,
      const arma::vec &i_a_t,
      const bool compute_H_and_z,
      const int bin_number,
      const double bin_tstart, const double bin_tstop);

public:
  EKF_solver(problem_data_EKF &p_, const std::string model_);

  void solve();
};

template<typename T>
class EKF_filter_worker{
  void do_comps(const arma::uvec::const_iterator it, int i,
                const arma::vec &i_a_t, const bool compute_z_and_H,
                const int bin_number,
                const double bin_tstart, const double bin_tstop);
  // Variables for computations
  problem_data_EKF &dat;
  arma::uvec::const_iterator first;
  const arma::uvec::const_iterator last;
  const arma::vec &i_a_t;
  const bool compute_z_and_H;
  const int i_start;
  const int bin_number;
  const double bin_tstart;
  const double bin_tstop;

  // local variables to compute temporary result
  arma::colvec u_;
  arma::mat U_;

public:
  EKF_filter_worker(
    problem_data_EKF &p_data,
    arma::uvec::const_iterator first_, const arma::uvec::const_iterator last_,
    const arma::vec &i_a_t_, const bool compute_z_and_H_,
    const int i_start_, const int bin_number_,
    const double bin_tstart_, const double bin_tstop_);

  void operator()();
};

struct EKF_filter_worker_calculations {
  /* (y_{it} - h(eta)) |_{\eta = x_{it}^T \alpha_{t}} */
  const double Y_residual;
  /* \frac{ \dot{h}(\eta) }{ Var(\eta)+ \xi } |_{\eta = x_{it}^T \alpha_{t}} */
  const double score_factor;
  /* \frac{ \dot{h}(\eta)^2 }{ Var(\eta)+ \xi } |_{\eta = x_{it}^T \alpha_{t}} */
  const double hessian_factor;
  /* \dot{h}(\eta) */
  const double var_inv;
  const double z_dot_factor;
};

struct EKF_logit_cals{
  static EKF_filter_worker_calculations cal(
      const bool do_die, const double time_outcome, const double at_risk_length,
      const double eta, const double denom_term);
};

struct EKF_exp_bin_cals {
  static EKF_filter_worker_calculations cal(
      const bool do_die, const double time_outcome, const double at_risk_length,
      const double eta, const double denom_term);
};

struct EKF_exp_clip_cals {
  static EKF_filter_worker_calculations  cal(
      const bool do_die, const double time_outcome, const double at_risk_length,
      const double eta, const double denom_term);
};

struct EKF_exp_clip_w_jump_cals{
  static EKF_filter_worker_calculations cal(
      const bool do_die, const double time_outcome, const double at_risk_length,
      const double eta, const double denom_term);
};




// UKF
class UKF_solver_Org : public Solver{
  problem_data &p_dat;
  const arma::uword m;
  const double k;
  const double w_0;
  const double w_i;
  const double sqrt_m_k;
  arma::mat sigma_points;

  inline void compute_sigma_points(const arma::vec &a_t,
                                   arma::mat &s_points,
                                   const arma::mat &P_x_x){
    arma::mat cholesky_decomp;
    if(!arma::chol(cholesky_decomp, P_x_x, "lower")){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Cholesky decomposition failed");
    }

    s_points.col(0) = a_t;
    for(arma::uword i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2);
  }

public:
  UKF_solver_Org(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &k_);

  void solve();
};





template<class T>
class UKF_solver_New : public Solver{
protected:
  problem_data &p_dat;
  const arma::uword m;
  const double a;
  const double k;
  const double b;
  const double lambda;
  const double w_0;
  const double w_0_c;
  const double w_0_cc;
  const double w_i;
  const double sqrt_m_lambda;
  arma::mat sigma_points;

  arma::vec weights_vec;
  arma::vec weights_vec_inv;
  arma::vec weights_vec_c;
  arma::vec weights_vec_c_inv;
  arma::vec weights_vec_cc;

  inline void compute_sigma_points(const arma::vec &a_t,
                                   arma::mat &s_points,
                                   const arma::mat &P_x_x){
    arma::mat cholesky_decomp;
    if(!arma::chol(cholesky_decomp, P_x_x, "lower")){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Cholesky decomposition failed");
    }

    s_points.col(0) = a_t;
    for(arma::uword i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_lambda * cholesky_decomp.unsafe_col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_lambda * cholesky_decomp.unsafe_col((i - 1) / 2);
  }

public:
  UKF_solver_New(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                 Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                 Rcpp::Nullable<Rcpp::NumericVector> &beta);

  void solve();
};

class UKF_solver_New_hepler_logit{
public:
  static constexpr bool need_risk_len = false;
  static constexpr bool adj_risk_len = false;

  static void mean_n_var_in_place(
      double&, const double, double &, const double);

  static double outcome(
      const bool, const double, const double, const double, const double);
};

class UKF_solver_New_hepler_exp_bin{
public:
  static constexpr bool need_risk_len = true;
  static constexpr bool adj_risk_len = false;

  static void mean_n_var_in_place(
      double&, const double, double &, const double);

  static double outcome(
      const bool, const double, const double, const double, const double);
};

class UKF_solver_New_hepler_exp_clip_time{
public:
  static constexpr bool need_risk_len = true;
  static constexpr bool adj_risk_len = true;

  static void mean_n_var_in_place(
      double&, const double, double &, const double);

  static double outcome(
      const bool, const double, const double, const double, const double);
};

class UKF_solver_New_hepler_exp_clip_time_w_jump{
public:
  static constexpr bool need_risk_len = true;
  static constexpr bool adj_risk_len = true;

  static void mean_n_var_in_place(
      double&, const double, double &, const double);

  static double outcome(
      const bool, const double, const double, const double, const double);
};

// Solver with approximation at the posterior mode sequentially

class SMA_hepler_logit{
private:
  static double NR_delta(
      const double offset, const double coef1, const double coef2,
      const double w, const double c0, const bool is_event,
      const double length);

public:
  static double compute_length(
      const double, const double, const double, const double,
      const bool, const double);

  static double second_d(const double, const double,
                         const double);
};

class SMA_hepler_exp{
private:
  static double NR_delta(
      const double offset, const double coef1, const double coef2,
      const double w, const double c0, const bool is_event,
      const double length);

public:
  static double compute_length(
      const double, const double, const double, const double, const bool,
      const double);

  static double second_d(const double, const double,
                         const double);
};

template<class T>
class SMA : public Solver
{
  problem_data &p_dat;
  std::string method;

public:
  SMA(problem_data &p_, std::string method_):
  p_dat(p_), method(method_)
  {
    if(method != "woodbury" && method != "cholesky")
      Rcpp::stop("Method '", method, "' not implemented");
  };

  void solve();
};

// Solver with approximation at the posterior mode sequentially

class GMA_hepler_logit{
public:
  static double d1(
      const double, const bool, const double);

  static double d2(
      const double, const double);
};

class GMA_hepler_exp{
public:
  static double d1(
      const double, const bool, const double);

  static double d2(
      const double, const double);
};

template<class T>
class GMA : public Solver
{
  problem_data &p_dat;
  const signed int max_rep;
  const double NR_eps;
  bool have_failed_once = false;

public:
  GMA(problem_data &p_, signed int max_rep_, double NR_eps_):
  p_dat(p_), max_rep(max_rep_), NR_eps(NR_eps_)
  { };

  void solve();
};

#endif
