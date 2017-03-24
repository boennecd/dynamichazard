#ifndef DDHAZARD_H
#define DDHAZARD_H

// [[Rcpp::plugins(cpp11)]]
#include "ddhazard_problem_data.h"
#include "arma_n_rcpp.h"



inline bool is_exponential_model(std::string model){
  return(model == "exp_combined" ||
         model == "exp_bin" ||
         model == "exp_clip_time" ||
         model == "exp_clip_time_w_jump");
}

// Abstact solver class
class Solver {
public:
  virtual void solve() = 0;
};

// Classes for EKF method
class EKF_filter_worker{
protected:
  virtual void do_comps(const arma::uvec::const_iterator it, int &i,
                        const arma::vec &i_a_t, const bool &compute_z_and_H,
                        const int &bin_number,
                        const double &bin_tstart, const double &bin_tstop) = 0; // abstact method to be implemented

  bool is_first_call;
  problem_data_EKF &dat;

  // local variables to compute temporary result
  arma::colvec u_;
  arma::mat U_;

public:
  EKF_filter_worker(problem_data_EKF &p_data);

  void operator()(arma::uvec::const_iterator first, const arma::uvec::const_iterator &last,
                const arma::vec &i_a_t, const bool &compute_z_and_H,
                const int &i_start, const int &bin_number,
                const double &bin_tstart, const double &bin_tstop);
};



class EKF_helper{
  unsigned long const max_threads;
  problem_data_EKF &p_data;
  std::vector<std::shared_ptr<EKF_filter_worker> > workers;
  const std::string model;

public:
  EKF_helper(problem_data_EKF &p_data_, const std::string model_);

  void parallel_filter_step(arma::uvec::const_iterator first, arma::uvec::const_iterator last,
                            const arma::vec &i_a_t,
                            const bool &compute_H_and_z,
                            const int &bin_number,
                            const double &bin_tstart, const double &bin_tstop);
};


class EKF_solver : public Solver{
  problem_data_EKF &p_dat;
  EKF_helper filter_helper;

public:
  EKF_solver(problem_data_EKF &p_, const std::string model);

  void solve();
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

  virtual void Compute_intermediates(const arma::uvec &r_set,
                                     const arma::vec offsets,
                                     const int t,
                                     const double bin_tstart, const double bin_tstop,
                                     arma::vec &c_vec, arma::mat &O) = 0;

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




class UKF_solver_New_logit : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O);

public:
  using UKF_solver_New::UKF_solver_New;
};


class UKF_solver_New_exp_bin : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O);

public:
  using UKF_solver_New::UKF_solver_New;
};


class UKF_solver_New_exp_clip_time : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O);

public:
  using UKF_solver_New::UKF_solver_New;
};


class UKF_solver_New_exp_clip_time_w_jump : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O);

public:
  using UKF_solver_New::UKF_solver_New;
};


class UKF_solver_New_exponential : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O);

public:
  using UKF_solver_New::UKF_solver_New;
};





// Solver with approximation at the posterior

class Posterior_approx_hepler_logit{
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




class Posterior_approx_hepler_exp{
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
class Posterior_approx : public Solver
{
  problem_data &p_dat;

public:
  Posterior_approx(problem_data &p_):
  p_dat(p_)
  {};

  void solve();
};

#endif
