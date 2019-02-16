#ifndef DDHAZARD_H
#define DDHAZARD_H

// [[Rcpp::plugins(cpp11)]]
#include "problem_data.h"
#include "arma_n_rcpp.h"
#include "family.h"

inline bool is_exponential_model(std::string model){
  return(model == "exp_bin" ||
         model == "exp_clip_time" ||
         model == "exp_clip_time_w_jump" ||
         model == "exponential");
}

// Abstact solver class
class Solver {
public:
  virtual void solve() = 0;

  // create a virtual, default destructor
  virtual ~Solver() = default;
};

// Classes for EKF method
class ddhazard_data_EKF;
class thread_pool;

class EKF_solver : public Solver {
  ddhazard_data &org;
  std::unique_ptr<ddhazard_data_EKF> p_dat;
  std::unique_ptr<thread_pool> pool;
  const std::string model;
  family_base &fam;

  void parallel_filter_step(
      arma::uvec::const_iterator first, arma::uvec::const_iterator last,
      const arma::vec &dynamic_coefs,
      const bool compute_H_and_z,
      const int bin_number,
      const double bin_tstart, const double bin_tstop);

public:
  EKF_solver(
    ddhazard_data&, const std::string,
    Rcpp::Nullable<Rcpp::NumericVector>, const unsigned int,
    const int, family_base &fam);

  void solve();
};

// UKF
class UKF_solver_Org : public Solver{
  ddhazard_data &p_dat;
  const arma::uword m;
  const double k;
  const double w_0;
  const double w_i;
  const double sqrt_m_k;
  arma::mat sigma_points;

  inline void compute_sigma_points(
      const arma::vec &a_t, arma::mat &s_points, const arma::mat &P_x_x){
    arma::mat cholesky_decomp;
    if(!arma::chol(cholesky_decomp, P_x_x, "lower")){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Cholesky decomposition failed");
    }

    s_points.col(0) = a_t;
    for(arma::uword i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_k * cholesky_decomp.col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_k * cholesky_decomp.col((i - 1) / 2);
  }

public:
  UKF_solver_Org(
    ddhazard_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &k_);

  void solve();
};





class UKF_solver_New : public Solver {
protected:
  ddhazard_data &p_dat;
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
  family_base &fam;

  void compute_sigma_points(const arma::vec&, arma::mat&, const arma::mat&);

public:
  UKF_solver_New(ddhazard_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                 Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                 Rcpp::Nullable<Rcpp::NumericVector> &beta,
                 family_base &fam);

  void solve();
};

// Solver with approximation at the posterior mode sequentially
class SMA : public Solver
{
  ddhazard_data &p_dat;
  std::string method;
  family_base &fam;

public:
  SMA(ddhazard_data &p_, std::string method_, family_base &fam):
  p_dat(p_), method(method_), fam(fam)
  {
    if(method != "woodbury" && method != "cholesky")
      Rcpp::stop("Method '", method, "' not implemented");
  };

  void solve();

  static double compute_length(
      const double, const double, const double,
      const double, const bool, const double, family_base &);
};

// Solver with approximation at the posterior mode sequentially
class GMA : public Solver
{
  ddhazard_data &p_dat;
  const unsigned int max_rep;
  const double NR_eps;
  bool have_failed_once = false;
  family_base &fam;

public:
  GMA(ddhazard_data &p, unsigned int max_rep, double NR_eps, family_base &fam):
  p_dat(p), max_rep((unsigned int)max_rep), NR_eps(NR_eps), fam(fam)
  { };

  void solve();
};


std::vector<double> logLike_cpp(
    const arma::mat&, const Rcpp::List&,
    const arma::mat&, const arma::mat&,
    arma::mat Q, const arma::mat&,
    const arma::vec&, const arma::vec&,
    const arma::vec&,
    const int, const std::string);

#endif
