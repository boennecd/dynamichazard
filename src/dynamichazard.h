// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <thread>
#include <future>
#include <algorithm>
#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#if defined(USE_OPEN_BLAS)
// Used to set the number of threads later
#include "cblas.h"
extern void openblas_set_num_threads(int num_threads);
extern int openblas_get_num_threads();
#define ARMA_USE_OPENBLAS
#define ARMA_DONT_USE_WRAPPER
#else
#define ARMA_USE_BLAS
#endif

// we know these are avialble with all R installations
#define ARMA_USE_LAPACK

#define ARMA_HAVE_STD_ISFINITE
#define ARMA_HAVE_STD_ISINF
#define ARMA_HAVE_STD_ISNAN
#define ARMA_HAVE_STD_SNPRINTF

// Rcpp has its own stream object which cooperates more nicely
// with R's i/o -- and as of Armadillo 2.4.3, we can use this
// stream object as well
#if !defined(ARMA_DEFAULT_OSTREAM)
#define ARMA_DEFAULT_OSTREAM Rcpp::Rcout
#endif

//#define ARMA_NO_DEBUG
// Note: This also disables the check in inv(A, B) for whether inversion is succesfull it seems http://arma.sourceforge.net/docs.html#inv
// from armadillo config.hpp
//// Uncomment the above line if you want to disable all run-time checks.
//// This will result in faster code, but you first need to make sure that your code runs correctly!
//// We strongly recommend to have the run-time checks enabled during development,
//// as this greatly aids in finding mistakes in your code, and hence speeds up development.
//// We recommend that run-time checks be disabled _only_ for the shipped version of your program.


#define ARMA_DONT_PRINT_ERRORS
// from armadillo config.hpp:
//// Comment out the above line if you don't want errors and warnings printed (eg. failed decompositions)

#include <RcppArmadillo.h> // has to come after defines: http://artax.karlin.mff.cuni.cz/r-help/library/RcppArmadillo/html/RcppArmadillo-package.html
// [[Rcpp::depends("RcppArmadillo")]]

// Classes for R like distributions families
class dist_family {
public:
  virtual double link_func(const double&) const = 0;
  virtual double link_func_inv(const double&) const = 0;
  virtual double variance(const double&) const = 0;
  virtual double d_mu_d_eta(const double&) const = 0; // d mu / d eta
  virtual double dev_resids(const double&, const double&, const double&) const = 0;

  virtual double time_offset(const double&) const = 0;
};

class logit_fam : public dist_family {
private:
  static constexpr double THRESH = 30.;
  static constexpr double MTHRESH = -30.;
  static constexpr double INVEPS = 1 / DOUBLE_EPS;

public:
  double link_func(const double&) const;
  double link_func_inv(const double&) const;
  double variance(const double&) const;
  double d_mu_d_eta(const double&) const; // d mu / d eta
  double dev_resids(const double&, const double&, const double&) const;

  double time_offset(const double&) const;
};

class poisson_fam : public dist_family
{
public:
  double link_func(const double&) const;
  double link_func_inv(const double&) const;
  double variance(const double&) const;
  double d_mu_d_eta(const double&) const; // d mu / d eta
  double dev_resids(const double&, const double&, const double&) const;

  double time_offset(const double&) const;
};


// Functions and classes for fixed effects. Similar to object in bigglm
int binomialCoeff(int n, int k);

class qr_obj{
public:
  qr_obj(unsigned int p):
  D(new arma::vec(p, arma::fill::zeros)), rbar(new arma::vec((p == 1)? 0 : binomialCoeff(p, 2), arma::fill::zeros)),
  thetab(new arma::vec(p, arma::fill::zeros)), ss(0.), checked(false),
  tol(new arma::vec(p, arma::fill::zeros))
  {}
  qr_obj() = default;

  std::shared_ptr<arma::vec> D;
  std::shared_ptr<arma::vec> rbar;
  std::shared_ptr<arma::vec> thetab;
  double ss;
  bool checked;
  std::shared_ptr<arma::vec> tol;
};


template<class T>
class bigglm_updateQR{
  // match logic in update.bigqr
  arma::vec linkinv(const arma::vec &eta);

  arma::vec d_mu_d_eta(const arma::vec &eta);

  arma::vec variance(const arma::vec &mu);

protected:
  T t;

public:
  bigglm_updateQR<T>(): t() {}

  void update(qr_obj &qr, // Previous/starting value. Will be overwritten
              const arma::mat &X, const arma::vec &eta,
              const arma::vec &offset, arma::vec &y, // y will not be altered
              const arma::vec &w);
};

using bigglm_updateQR_logit   = bigglm_updateQR<logit_fam>;
using bigglm_updateQR_poisson = bigglm_updateQR<poisson_fam>;

arma::vec bigglm_regcf(qr_obj &qr);

bool is_exponential_model(std::string model);






// problem data classes
// Hepler structure to reference data

#ifndef DDHAZARD_PROBLEM_DATA
#define DDHAZARD_PROBLEM_DATA

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

  const double ridge_eps;

  const bool debug;
  const double LR;

  const bool any_fixed_terms_in_state_vec;
  const bool use_pinv;
  const std::string criteria;

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
               const double ridge_eps_,
               const bool use_pinv_,
               const std::string criteria_):
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
    ridge_eps(ridge_eps_),

    debug(debug_),
    LR(LR_.isNotNull() ? Rcpp::as< Rcpp::NumericVector >(LR_)[0] : 1.0),

    any_fixed_terms_in_state_vec(n_fixed_terms_in_state_vec_ > 0),
    use_pinv(use_pinv_),
    criteria(criteria_),

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
                   const double ridge_eps_,
                   const bool use_pinv_,
                   const std::string criteria_):
    problem_data(n_fixed_terms_in_state_vec_, X, fixed_terms, tstart_, tstop_, is_event_in_bin_, a_0,
                 fixed_parems_start, Q_0_, Q_, risk_obj, F__,
                 eps_fixed_parems_, max_it_fixed_params_,
                 weights_,
                 n_max, eps, verbose,
                 order_, est_Q_0, debug_,
                 LR_,
                 n_threads_, ridge_eps_,
                 use_pinv_, criteria_),

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






#ifndef DD_UTIL_FUNCS
#define DD_UTIL_FUNCS

// Print function to print out vectors as rows
template<typename T>
inline
  void
  my_print(const T& X, std::string msg = "")
  {
    if(msg != "")
      Rcpp::Rcout << msg << std::endl;

    if(X.n_cols > 1){
      X.print();

    } else{
      for(arma::uword col=0; col < X.n_cols; ++col)
      {
        for(arma::uword row=0; row < X.n_rows; ++row)
          Rcpp::Rcout << std::setw(10) << std::setprecision(5) <<  X(row,col) << ' ';

        Rcpp::Rcout << std::endl;
      }
    }
  }



// Inversion methods
template<typename eT, typename T2>
inline
  void inv(arma::Mat<eT>& out, T2&& X,
           const bool use_pinv = false,
           const std::string err_msg = ""){
    if(use_pinv){
      // Compute the pseudo inverse using SVD
      // Tolerance controls the value for which eigenvalues and vectors are
      // dropped
      // the default tolerance is max(m,n)*max_sv*datum::eps, where:
      //   m = number of rows and n = number of columns
      //   max_sv = maximal singular value
      //   datum::eps = difference between 1 and the least value greater than 1 that is representable

      // Two method are avialable for the SVD: "dc" and "std" (former is default)
      // See armadillo-#.###.#\include\armadillo_bits\op_pinv_meat.hpp
      // "dc" uses auxlib::svd_dc_econ which calls lapack::cx_gesdd
      // "std" uses auxlib::svd_econ which calls lapack::cx_gesvd
      // see armadillo-#.###.#\include\armadillo_bits\auxlib_meat.hpp
      // and armadillo-#.###.#\include\armadillo_bits\wrapper_lapack.hpp

      // For double, former calls zgesdd and latter calls arma_zgesvd

      // Acording to this post: https://groups.google.com/forum/#!topic/julia-dev/mmgO65i6-fA
      // "The routine dgesdd uses a divide and conquer algorithm for the
      //   bidiagonal SVD, whereas dgesvd uses a QR algorithm. The former is
      //   faster in cases where there is significant deflation and should be
      //   generally preferred for large matrices."

      // A next question is regarding whether the tolerance. Wikipedia suggest
      // this tolerance is the default in other langagues as well https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse
      if(!arma::pinv(out, std::forward<T2>(X))){
        Rcpp::stop(err_msg);
      }
    } else{
      if(!arma::inv(out, std::forward<T2>(X))){
        Rcpp::stop(err_msg);
      }
    }
  }

template<typename eT, typename T2>
inline
  void inv_sympd(arma::Mat<eT>& out, T2&& X,
                 const bool use_pinv = false,
                 const std::string err_msg = ""){

    if(use_pinv){
      inv(out, std::forward<T2>(X), true, err_msg);
    } else if(!arma::inv_sympd(out, std::forward<T2>(X))){
      Rcpp::stop(err_msg);
    }
  }


#endif
