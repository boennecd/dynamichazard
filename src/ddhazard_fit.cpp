#include "dynamichazard.h"
#include "thread_pool.h"

using uword = arma::uword;

bool is_exponential_model(std::string model){
  return(model == "exp_combined" ||
         model == "exp_bin" ||
         model == "exp_trunc_time" ||
         model == "exp_trunc_time_w_jump");
}

// Define convergence criteria
inline double relative_norm_change(const arma::mat &prev_est, const arma::mat &new_est){
  return arma::norm(prev_est - new_est, 2) / (arma::norm(prev_est, 2) + 1.0e-10);
}
double (*conv_criteria)(const arma::mat&, const arma::mat&) = relative_norm_change;

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

// Hepler structure to reference data
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
               const int n_max = 100, const double eps = 0.001,
               const bool verbose = false,
               const int order_ = 1, const bool est_Q_0 = true,
               const bool debug_ = false,
               Rcpp::Nullable<Rcpp::NumericVector> LR_ = R_NilValue,
               const int n_threads_ = -1,
               const double ridge_eps_ = .0001):
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
                   const int n_max = 100, const double eps = 0.001,
                   const bool verbose = false,
                   const int order_ = 1, const bool est_Q_0 = true,
                   const bool is_cont_time = false,
                   const unsigned int NR_it_max_ = 1000,
                   const bool debug_ = false,
                   const int n_threads_ = -1,
                   const double ridge_eps_ = .0001):
    problem_data(n_fixed_terms_in_state_vec_, X, fixed_terms, tstart_, tstop_, is_event_in_bin_, a_0,
                 fixed_parems_start, Q_0_, Q_, risk_obj, F__,
                 eps_fixed_parems_, max_it_fixed_params_,
                 n_max, eps, verbose,
                 order_, est_Q_0, debug_,
                 LR_,
                 n_threads_, ridge_eps_),

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

// Abstact solver class
class Solver {
public:
  virtual void solve() = 0;
};





namespace exp_model_funcs {
// Namespace to avoid namespace pollution and avoid error-40 for Taylor/Laruens series
// By definition:
//  eps                 epsilon in ridge regression like solution
//  a                   at risk length
//  eta                 linear predictor x^T * beta
//  v                   a * exp(eta)
//  exp_x               exp(x)
//  inv_x               inv(x)
//  expect_chance_die   1 - exp(- v)

inline double expect_time(const double v, const double a,
                          const double inv_exp_v, const double exp_eta){
  return((v >= 1e-8) ? (1.0 - inv_exp_v) / exp_eta :
           a + exp_eta*(-pow(a,2)/2 + exp_eta*(pow(a,3)/6 - (pow(a,4)*exp_eta)/24))
  );
}

inline double expect_time_w_jump(const double exp_eta, const double inv_exp_eta,
                                 const double inv_exp_v, const double a){
  return((a * exp_eta >= 1e-4) ?
           (inv_exp_v - 1)*(-1 + a*exp_eta)*inv_exp_eta :
           a - (3*pow(a,2)*exp_eta)/2);
}

inline double expect_chance_die(const double v, const double inv_exp_v){
  return((v >= 1e-5) ? 1.0 - inv_exp_v :
           v * (1.0 - v / 2.0 * (1.0 + v / 6.0 * (1.0 - v / 24 * (1.0 - v /120.0)))));
}

inline double inv_var_wait_time(const double v, const double exp_eta, const double inv_exp_v){
  return((v >= 1e-4) ?
           // Set v = delta * exp(eta)
           // Then: exp(2eta) / (1 - exp(-2 delta * exp(eta)) - 2 exp(-delta * exp(eta)) delta exp(eta)) =
           //               exp(2eta) / (1 - exp(-2v) - 2 * v * exp(-v))
           exp_eta * exp_eta / (1.0 - inv_exp_v * inv_exp_v - 2.0 * v * inv_exp_v) :
           // Laruent series from https://www.wolframalpha.com/input/?i=1%2F(1-exp(2v)-2v*exp(v))
           exp_eta * exp_eta *
             (-1 / v * (1 / 4 - v * (1 / 4 - v * (5 / 48 - v * (1/48 - v /1440)))))
  );
}

inline double inv_var_chance_die(const double v, const double inv_exp_v){
  return((v >= 1e-4) ?
           // Set v = a exp(eta)
           // Then: 1 / exp(- a exp(eta))(1 - exp(-a exp(eta))) = 1 / exp(-v) (1 - exp(-v))
           1 / (inv_exp_v * (1 - inv_exp_v)) :
           //Lauren series from https://www.wolframalpha.com/input/?i=1%2F((1-exp(-v))exp(-v))
           1 / v * (1 + v * (3 / 2 + v * (13 / 12 + v * (1 / 2 + v * 119 / 720))))
  );

}

inline double EKF_fac_score_die(const double exp_eta, const double v,
                                const double exp_v, const double a,
                                const double eps){
  if(a * exp_eta >= 1e-4){
    return(-(((-1 + exp_eta + exp_v*(-1 + a*(-2 + exp_eta))*exp_eta +
           pow(exp_v,2)*(1 + eps*pow(exp_eta,2)))*v)/
             (-1 + pow(exp_v,2)*(-1 + 2*eps*v - eps*pow(exp_eta,2)) -
               pow(exp_v,3)*eps*(1 + eps*pow(exp_eta,2)) + exp_v*(2 + eps + pow(v,2) +
               eps*pow(exp_eta,2)))));
  } else {
    // Taylor expansion around exp_eta = 0
    const double v_eps_ratio = v / eps;

    return(v_eps_ratio * (
        1 + v_eps_ratio * (
            (-2 + a - 2*eps) / 2 - v_eps_ratio *
              (-12 + 6*a - 3*pow(a,2) + 2*pow(a,3) - 30*eps + 14*a*eps - 6*pow(eps,2)) / 12)));
  }
}

inline double EKF_fac_score_time(const double exp_eta, const double v,
                                 const double exp_v, const double a,
                                 const double eps){
  if(a * exp_eta >= 1e-3){
    return(-(((-1 + exp_v - v)*(1 - exp_eta + eps*exp_eta*pow(exp_v,2) + exp_v*(-1 + exp_eta + v)))/
           (1 - (2 + eps + pow(v,2) + eps*pow(exp_eta,2))*exp_v + eps*(1 + eps*pow(exp_eta,2))*pow(exp_v,3) +
             pow(exp_v,2)*(1 + eps*pow(exp_eta,2) - 2*eps*v))));
  } else{
    const double exp_eta_eps_ratio = exp_eta / eps;

    return(exp_eta_eps_ratio * (-pow(a,2)/2 + exp_eta_eps_ratio *
               (-3*pow(a,4) + 2*pow(a,5) + 4*pow(a,3)*eps)/12));
  }
}

inline double EKF_fac_var(const double exp_eta, const double v,
                          const double exp_v, const double a,
                          const double eps){
  if(a * exp_eta >= 1e-3){
    return((-1 + exp_v*(3 + pow(v,2) + eps*pow(exp_v,3) + exp_v*(-3 + eps + pow(v,2)*(-1 + eps) +
           pow(v,2)*eps*pow(exp_eta,2)+ 2*eps*v) + pow(exp_v,2)*(1 - 2*eps*(1 + v))))/
             (exp_v - (2 + eps + (pow(a,2) + eps)*pow(exp_eta,2))*pow(exp_v,2) + (1 + eps*exp_eta*(-2*a + exp_eta))*
               pow(exp_v,3) + eps*(1 + eps*pow(exp_eta,2))*pow(exp_v,4)));
  } else{
    return((4*pow(a,2) + pow(a,4))*exp_eta * (exp_eta / (4 * eps)));
  }
}

inline double var_wait_time(const double v, const double a,  const double exp_eta, const double inv_exp_v){
  // exp(eta)^(-2) * (1 - exp(-2v) - 2 * v * exp(-v))
  // v = a * exp(eta) => ... = a^2 * v^(-2) * (1 - exp(-2v) - 2 * v * exp(-v))
  //
  // Use Taylor series for the latter when v is small: https://www.wolframalpha.com/input/?i=(1-exp(-2v)-2v*exp(-v))%2F(v%5E2)

  return((v >= 1e-4) ?
           (1.0 - inv_exp_v * inv_exp_v - 2.0 * v * inv_exp_v) / (exp_eta * exp_eta) :
           a * a * v * (1/3 - v * (1/3 - v * (11/60 - v * (13/180 - v * 19/840)))));
}

inline double var_wait_time_w_jump(const double exp_eta,
                                   const double inv_exp_v, const double a){
  return((a * exp_eta >= 1e-4) ?
           (1 - pow(inv_exp_v,2) + pow(a*exp_eta,2)*(3*inv_exp_v - pow(inv_exp_v,2)) +
              a*exp_eta*(-4*inv_exp_v + 2*pow(inv_exp_v,2)))/ pow(exp_eta,2) :
           exp_eta*((7*pow(a,3))/3 + exp_eta*((-19*pow(a,4))/6 + exp_eta*(
               (34*pow(a,5))/15 - (407*pow(a,6)*exp_eta)/360)))
  );
}

inline double var_chance_die(const double v, const double inv_exp_v){
  // Taylor series from https://www.wolframalpha.com/input/?i=exp(-v)+*+(1+-+exp(-v))
  return((v >= 1e-4) ?
           inv_exp_v * (1 - inv_exp_v) :
           v * (1 - v * (3/2 - v * (7/6 - v * (5/8 - v * 31 /120)))));
}

inline double covar(const double v,
                    const double inv_exp_v, const double exp_eta){
  // Taylor series from http://www.wolframalpha.com/input/?i=-1+*+exp(-2v)+*+(1+%2B+v+*+exp(v)+-+exp(v))
  return(
    (v >= 1e-4) ?
  -1 * inv_exp_v * (inv_exp_v + v - 1.) / exp_eta :
  -v * v * (1/2 - v * (2/3 - v * (11/24 - v * (13/60 - 19 * v / 240)))) / exp_eta);
}

inline double binary_score_fac(const double v, const double  inv_exp_v, double eps){
  return((inv_exp_v*v)/(eps + (1 - inv_exp_v)*inv_exp_v));
}

inline double binary_var_fac(const double v, const double  inv_exp_v, double eps){
  return(pow(inv_exp_v*v, 2)/(eps + (1 - inv_exp_v)*inv_exp_v));
}


inline double trunc_time_score_fac(const double exp_eta, const double inv_exp_eta,
                                        const double inv_exp_v, const double a,
                                        const double eps){

  return((a * exp_eta >= 1e-4) ?
           (-1 + inv_exp_v + a*exp_eta*inv_exp_v)/(eps*exp_eta + inv_exp_eta - 2*a*inv_exp_v -
                                                    inv_exp_eta*pow(inv_exp_v,2)) :
           ((exp_eta / eps) *(-3*pow(a,2)*eps + (pow(a,5) + 2*pow(a,3)*eps)*exp_eta))/(6*eps));
}

inline double trunc_time_var_fac(const double exp_eta, const double inv_exp_eta,
                                      const double inv_exp_v, const double a,
                                      const double eps){
  return((a * exp_eta >= 1e-4) ?
           (1 - 2*inv_exp_v - 2*a*exp_eta*inv_exp_v + pow(inv_exp_v,2) + 2*a*exp_eta*pow(inv_exp_v,2) +
              pow(a*exp_eta*inv_exp_v, 2))/
                (1 + eps*pow(exp_eta,2) - 2*a*exp_eta*inv_exp_v - pow(inv_exp_v,2)) :
           (pow(exp_eta/eps,2)*(3*pow(a,4)*eps + (-pow(a,7) - 4*pow(a,5)*eps)*exp_eta))/12);
}

inline double trunc_time_w_jump_score_fac(const double exp_eta, const double v,
                                               const double inv_exp_v, const double a,
                                               const double eps){
  if(a * exp_eta >= 1e-5){
    return((exp_eta*(-1 + inv_exp_v) + a*pow(exp_eta,2)*inv_exp_v - pow(a*exp_eta,2) *exp_eta*inv_exp_v)/
           (1 - pow(inv_exp_v,2) + a*exp_eta*(-4*inv_exp_v + 2*pow(inv_exp_v,2)) +
             pow(exp_eta,2)*(eps + pow(a,2)*(3*inv_exp_v - pow(inv_exp_v,2)))));
  } else{
    return((pow(v,2)*(-18*pow(a,2) - 9*eps + (7*pow(a,2) + 8*eps)*v))/(24*pow(a,4)*exp_eta + 24*pow(a,2)*eps*exp_eta +
           6*pow(eps,2)*exp_eta));
  }
}

inline double trunc_time_w_jump_var_fac(const double exp_eta,
                                             const double inv_exp_v, const double a,
                                             const double eps){
  if(a * exp_eta >= 1e-4){
    return((1 - 2*inv_exp_v + pow(inv_exp_v,2) + (- 2*pow(a*exp_eta,3) + pow(a*exp_eta,4))*pow(inv_exp_v,2) +
           pow(a*exp_eta,2)*(2*inv_exp_v - pow(inv_exp_v,2)) + a*exp_eta*(-2*inv_exp_v + 2*pow(inv_exp_v,2)))/
             (1 - pow(inv_exp_v,2) + a*exp_eta*(-4*inv_exp_v + 2*pow(inv_exp_v,2)) +
               pow(exp_eta,2)*(eps + pow(a,2)*(3*inv_exp_v - pow(inv_exp_v,2)))));
  } else{
    return((9*pow(a,4)*pow(exp_eta,2))/(4*eps));
  }
}}

class EKF_helper{
  // worker class for parallel computation
  // This class is abstact as the method do_computation will differ between
  // the models
  class filter_worker{
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
    filter_worker(problem_data_EKF &p_data):
    is_first_call(true), dat(p_data)
    {}

    void operator()(arma::uvec::const_iterator first, const arma::uvec::const_iterator &last,
                  const arma::vec &i_a_t, const bool &compute_z_and_H,
                  const int &i_start, const int &bin_number,
                  const double &bin_tstart, const double &bin_tstop){
      // potentially intialize variables and set entries to zeroes in any case
      if(is_first_call){
        u_ = arma::vec(dat.n_params_state_vec);
        U_ = arma::mat(dat.n_params_state_vec, dat.n_params_state_vec);
        is_first_call = false;
      }
      u_.zeros();
      U_.zeros();

      // compute local results
      int i = i_start;
      for(arma::uvec::const_iterator it = first; it != last; it++){
        do_comps(it, i, i_a_t, compute_z_and_H, bin_number,
                 bin_tstart, bin_tstop);
      }

      // Update shared variable
      {
        std::lock_guard<std::mutex> lk(dat.m_U);
        dat.U(dat.span_current_cov, dat.span_current_cov) +=  U_;
      }

      {
        std::lock_guard<std::mutex> lk(dat.m_u);
        dat.u(dat.span_current_cov) += u_;
      }
    }
  };

  // worker for the logit model
  class filter_worker_logit : public filter_worker {
  private:
    void do_comps(const arma::uvec::const_iterator it, int &i,
                  const arma::vec &i_a_t, const bool &compute_z_and_H,
                  const int &bin_number,
                  const double &bin_tstart, const double &bin_tstop){
      const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);
      double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
      const double exp_eta = exp(arma::dot(i_a_t, x_) + offset);

      // Can be issue here with overflow in denominator
      double tmp_denom = pow(1.0 + exp_eta, 2.0);
      const double var = std::isinf(tmp_denom) ?
        pow(exp_eta, -1) : (exp_eta / pow(exp_eta + 1.0, 2.0));

      double score_fac = exp_eta/(dat.ridge_eps + exp_eta +
                                  2*dat.ridge_eps*exp_eta + dat.ridge_eps*pow(exp_eta,2));
      double var_fac = pow(exp_eta,2)/(pow(1 + exp_eta,4)*(dat.ridge_eps + exp_eta/pow(1 + exp_eta, 2)));

      u_ += x_ * (score_fac *
        ((dat.is_event_in_bin(*it) == bin_number) - exp_eta / (1.0 + exp_eta)));
      U_ += x_ *  (x_.t() * var_fac);

      if(compute_z_and_H){
        dat.H_diag_inv(i) = pow(var, -1);
        dat.z_dot(dat.span_current_cov, i) = x_ *  var;
        ++i;
      }
    }

  public:
    filter_worker_logit(problem_data_EKF &p_data):
    filter_worker(p_data)
    {}
  };

  // worker for the continous model with exponential distribution
  class filter_worker_exponential : public filter_worker {
  private:
    void do_comps(const arma::uvec::const_iterator it, int &i,
                  const arma::vec &i_a_t, const bool &compute_z_and_H,
                  const int &bin_number,
                  const double &bin_tstart, const double &bin_tstop){
      // Compute intermediates
      const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);

      double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
      const double eta = arma::dot(i_a_t, x_) + offset;

      const double do_die = (dat.is_event_in_bin(*it) == bin_number);
      const double time_outcome = std::min(dat.tstop(*it), bin_tstop) - std::max(dat.tstart(*it), bin_tstart);
      const double at_risk_length = do_die ? bin_tstop - std::max(dat.tstart(*it), bin_tstart) : time_outcome;

      const double exp_eta = exp(eta);
      const double inv_exp_eta = pow(exp_eta, -1);

      const double v = at_risk_length * exp_eta;
      const double exp_v = exp(v);
      const double inv_exp_v = pow(exp_v, -1.0);

      const double expect_time = exp_model_funcs::expect_time(
        v, at_risk_length, inv_exp_v, exp_eta);

      const double expect_chance_die = exp_model_funcs::expect_chance_die(v, inv_exp_v);

      const double fac_score_die = exp_model_funcs::EKF_fac_score_die(
        exp_eta, v, exp_v, at_risk_length, dat.ridge_eps);
      const double fac_score_time = exp_model_funcs::EKF_fac_score_time(
        exp_eta, v, exp_v, at_risk_length, dat.ridge_eps);
      const double var_fac = exp_model_funcs::EKF_fac_var(
        exp_eta, v, exp_v, at_risk_length, dat.ridge_eps);

      u_ += x_ * (
        fac_score_time * (time_outcome - expect_time) +
        fac_score_die * (do_die - expect_chance_die));

      U_ += x_ * (x_.t() * var_fac);

      if(compute_z_and_H){
        // Compute terms from waiting time
        dat.H_diag_inv(i) = exp_model_funcs::inv_var_wait_time(v, exp_eta, inv_exp_v);

        dat.z_dot(dat.span_current_cov, i) =  x_ * ((v >= 1e-6) ?
                                                inv_exp_v * (inv_exp_eta + at_risk_length) - inv_exp_eta :
                                                // Taylor series from https://www.wolframalpha.com/input/?i=exp(-v)%2Bv*exp(-v)-1
                                                inv_exp_eta * (- v * v) * (1/2 - v * (1/3 - v * (1/8 - v * (1/30 - v/144)))));

        // Compute terms from binary out come
        dat.H_diag_inv(i + dat.n_in_last_set) = exp_model_funcs::inv_var_chance_die(
          v, expect_chance_die);

        dat.z_dot(dat.span_current_cov, i + dat.n_in_last_set) =
          x_ * (at_risk_length * exp_eta * inv_exp_v);
        ++i;
      }
    }

  public:
    filter_worker_exponential(problem_data_EKF &p_data):
    filter_worker(p_data)
      {}
  };

  // worker for the continous model with exponential distribution where only the
  // binary variable is used
  class filter_worker_exp_bin : public filter_worker {
  private:
    void do_comps(const arma::uvec::const_iterator it, int &i,
                  const arma::vec &i_a_t, const bool &compute_z_and_H,
                  const int &bin_number,
                  const double &bin_tstart, const double &bin_tstop){
      // Compute intermediates
      const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);

      double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
      const double eta = arma::dot(i_a_t, x_) + offset;

      const double do_die = (dat.is_event_in_bin(*it) == bin_number);
      const double at_risk_length = std::min(dat.tstop(*it), bin_tstop) - std::max(dat.tstart(*it), bin_tstart);

      const double exp_eta = exp(eta);
      const double v = at_risk_length * exp_eta;
      const double exp_v = exp(v);
      const double inv_exp_v = pow(exp_v, -1.0);

      const double expect_chance_die = exp_model_funcs::expect_chance_die(v, inv_exp_v);

      const double score_fac = exp_model_funcs::binary_score_fac(v, inv_exp_v, dat.ridge_eps);
      const double var_fac = exp_model_funcs::binary_var_fac(v, inv_exp_v, dat.ridge_eps);

      u_ += x_ * (score_fac * (do_die - expect_chance_die));

      U_ += x_ * (x_.t() * var_fac);

      if(compute_z_and_H){
        // Compute terms from waiting time
        dat.H_diag_inv(i) = exp_model_funcs::inv_var_chance_die(v, expect_chance_die);

        dat.z_dot(dat.span_current_cov, i) =  x_ * (inv_exp_v * v);
        ++i;
      }
    }
  public:
    filter_worker_exp_bin(problem_data_EKF &p_data):
    filter_worker(p_data)
    {}
  };

  // worker for the continous model with exponential distribution where only the
  // right truncated variable is used
  class filter_worker_exp_trunc_time : public filter_worker {
  private:
    void do_comps(const arma::uvec::const_iterator it, int &i,
                  const arma::vec &i_a_t, const bool &compute_z_and_H,
                  const int &bin_number,
                  const double &bin_tstart, const double &bin_tstop){
      // Compute intermediates
      const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);

      double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
      const double eta = arma::dot(i_a_t, x_) + offset;

      const double do_die = (dat.is_event_in_bin(*it) == bin_number);
      const double time_outcome = std::min(dat.tstop(*it), bin_tstop) - std::max(dat.tstart(*it), bin_tstart);
      const double at_risk_length = do_die ? bin_tstop - std::max(dat.tstart(*it), bin_tstart) : time_outcome;

      const double exp_eta = exp(eta);
      const double inv_exp_eta = pow(exp_eta, -1);
      const double v = at_risk_length * exp_eta;
      const double exp_v = exp(v);
      const double inv_exp_v = pow(exp_v, -1.0);

      const double expect_time = exp_model_funcs::expect_time(v, at_risk_length, inv_exp_v, exp_eta);

      const double score_fac = exp_model_funcs::trunc_time_score_fac(
        exp_eta, inv_exp_eta, inv_exp_v, at_risk_length, dat.ridge_eps);
      const double var_fac = exp_model_funcs::trunc_time_var_fac(
        exp_eta, inv_exp_eta, inv_exp_v, at_risk_length, dat.ridge_eps);


      u_ += x_ * (score_fac * (time_outcome - expect_time));

      U_ += x_ * (x_.t() * var_fac);

      if(compute_z_and_H){
        // Compute terms from waiting time
        dat.H_diag_inv(i) = exp_model_funcs::inv_var_wait_time(v, exp_eta, inv_exp_v);

        dat.z_dot(dat.span_current_cov, i) =  x_ * (inv_exp_eta*inv_exp_v*(1 - exp_v + v));
        ++i;
      }
    }
  public:
    filter_worker_exp_trunc_time(problem_data_EKF &p_data):
    filter_worker(p_data)
    {}
  };

  // worker for the continous model with exponential distribution where only the
  // right truncated variable is used where outcomes are negative
  class filter_worker_exp_trunc_time_w_jump : public filter_worker {
  private:
    void do_comps(const arma::uvec::const_iterator it, int &i,
                  const arma::vec &i_a_t, const bool &compute_z_and_H,
                  const int &bin_number,
                  const double &bin_tstart, const double &bin_tstop){
      // Compute intermediates
      const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);

      double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
      const double eta = arma::dot(i_a_t, x_) + offset;

      const double do_die = (dat.is_event_in_bin(*it) == bin_number);
      const double time_outcome = std::min(dat.tstop(*it), bin_tstop) - std::max(dat.tstart(*it), bin_tstart);
      const double at_risk_length = do_die ? bin_tstop - std::max(dat.tstart(*it), bin_tstart) : time_outcome;

      const double exp_eta = exp(eta);
      const double inv_exp_eta = pow(exp_eta, -1);
      const double v = at_risk_length * exp_eta;
      const double exp_v = exp(v);
      const double inv_exp_v = pow(exp_v, -1);

      const double expect_time = exp_model_funcs::expect_time_w_jump(exp_eta, inv_exp_eta, inv_exp_v, at_risk_length);

      const double score_fac = exp_model_funcs::trunc_time_w_jump_score_fac(
        exp_eta, v, inv_exp_v, at_risk_length, dat.ridge_eps);
      const double var_fac = exp_model_funcs::trunc_time_w_jump_var_fac(
        exp_eta, inv_exp_v, at_risk_length, dat.ridge_eps);

      if(abs(score_fac) > 1e1 || var_fac <= 0){
        {
          std::lock_guard<std::mutex> lk(dat.m_u);
          Rcpp::Rcout << "exp_eta " << exp_eta
                      << "\tscore_fac " << score_fac
                      << "\tat_risk_length " << at_risk_length
                      << "\tvar_fac " << var_fac << std::endl;
        }
      }

      if(do_die){
        u_ += x_ * (score_fac * (time_outcome - expect_time - at_risk_length));
      } else {
        u_ += x_ * (score_fac * (time_outcome - expect_time));
      }

      U_ += x_ * (x_.t() * var_fac);

      if(compute_z_and_H){
        // Compute terms from waiting time
        dat.H_diag_inv(i) =
          1 / exp_model_funcs::var_wait_time_w_jump(exp_eta, inv_exp_v, inv_exp_v);

        dat.z_dot(dat.span_current_cov, i) =  x_ * (
          -inv_exp_eta + at_risk_length*inv_exp_v - pow(at_risk_length,2)*exp_eta*inv_exp_v + inv_exp_eta*inv_exp_v
        );
        ++i;
      }
    }
  public:
    filter_worker_exp_trunc_time_w_jump(problem_data_EKF &p_data):
    filter_worker(p_data)
    {}
  };

  unsigned long const max_threads;
  problem_data_EKF &p_data;
  std::vector<std::shared_ptr<filter_worker> > workers;
  const std::string model;

public:
  EKF_helper(problem_data_EKF &p_data_, const std::string model_):
  max_threads((p_data_.n_threads > 1) ? p_data_.n_threads - 1 : 1),
  p_data(p_data_), workers(), model(model_)
  {
    if(p_data.debug)
      Rcpp::Rcout << "EKF solver will use at most " << max_threads << " threads" << std::endl;
  }

  void parallel_filter_step(arma::uvec::const_iterator first, arma::uvec::const_iterator last,
                            const arma::vec &i_a_t,
                            const bool &compute_H_and_z,
                            const int &bin_number,
                            const double &bin_tstart, const double &bin_tstop){
    // Set entries to zero
    p_data.U.zeros();
    p_data.u.zeros();
    if(compute_H_and_z){
      p_data.z_dot.zeros();
      p_data.H_diag_inv.zeros();
    }

    // Compute the number of blocks to create
    unsigned long const length = std::distance(first, last);

    unsigned long const block_size = 250;
    unsigned long const num_blocks=(length+block_size-1)/block_size;
    std::vector<std::future<void> > futures(num_blocks-1);
    thread_pool pool(num_blocks - 1, max_threads);

    // Create workers if needed
    for(auto i = workers.size(); i < num_blocks; i++){
      if(model == "logit"){
        std::shared_ptr<filter_worker> new_p(new filter_worker_logit(p_data));
        workers.push_back(std::move(new_p));

      } else if (model == "exp_combined"){
        std::shared_ptr<filter_worker> new_p(new filter_worker_exponential(p_data));
        workers.push_back(std::move(new_p));

      } else if (model == "exp_bin"){
        std::shared_ptr<filter_worker> new_p(new filter_worker_exp_bin(p_data));
        workers.push_back(std::move(new_p));
      } else if(model == "exp_trunc_time"){
        std::shared_ptr<filter_worker> new_p(new filter_worker_exp_trunc_time(p_data));
        workers.push_back(std::move(new_p));
      } else if(model == "exp_trunc_time_w_jump"){
        std::shared_ptr<filter_worker> new_p(new filter_worker_exp_trunc_time_w_jump(p_data));
        workers.push_back(std::move(new_p));
      } else
        Rcpp::stop("EKF is not implemented for model '" + model  +"'");
    }

    // start workers
    // declare outsite of loop to ref after loop
    arma::uvec::const_iterator block_start = first;
    auto it = workers.begin();
    int i_start = 0;

    for(unsigned long i = 0; i < num_blocks - 1; ++i, ++it)
    {
      arma::uvec::const_iterator block_end = block_start;
      std::advance(block_end, block_size);

      auto func =
        [it, block_start, block_end, &i_a_t, &compute_H_and_z, i_start, &bin_number, &bin_tstart, &bin_tstop](){
          (*it->get())(block_start, block_end, i_a_t, compute_H_and_z,
           i_start, bin_number, bin_tstart, bin_tstop);
        };

        futures[i] = pool.submit(func);
        i_start += block_size;
        block_start = block_end;
    }
    (*(it->get()))(block_start, last, i_a_t, compute_H_and_z, i_start, bin_number, bin_tstart, bin_tstop); // compute last enteries on this thread

    for(unsigned long i = 0; i < num_blocks - 1; ++i)
    {
      futures[i].get();   // will throw if any of the threads did
    }
  }
};

class EKF_solver : public Solver{
  problem_data_EKF &p_dat;
  EKF_helper filter_helper;


public:
  EKF_solver(problem_data_EKF &p_, const std::string model):
  p_dat(p_), filter_helper(p_, model)
  {}

  void solve(){
    double bin_tstop = p_dat.min_start;
    for (int t = 1; t < p_dat.d + 1; t++){

      double bin_tstart = bin_tstop;
      double delta_t = p_dat.I_len[t - 1];
      bin_tstop += delta_t;

      // E-step: Filter step
      p_dat.a_t_less_s.col(t - 1) = p_dat.F_ *  p_dat.a_t_t_s.unsafe_col(t - 1);
      p_dat.V_t_less_s.slice(t - 1) = p_dat.F_ * p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ + delta_t * p_dat.Q;

      // E-step: scoring step: information matrix and scoring vector
      arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;
      arma::vec i_a_t = p_dat.a_t_less_s.col(t - 1);
      arma::mat V_t_less_s_inv;
      unsigned int n_NR_it = 0;

      while(true){
        ++n_NR_it;

#if defined(USE_OPEN_BLAS)
        openblas_set_num_threads(1);
#endif
        filter_helper.parallel_filter_step(r_set.begin(), r_set.end(), i_a_t(p_dat.span_current_cov), t == p_dat.d, t - 1,
                                           bin_tstart, bin_tstop);

        if(p_dat.u.has_inf() || p_dat.u.has_nan()){
          Rcpp::stop("ddhazard_fit_cpp estimation error: Score vector in correction step had inf or nan elements in bin " +
            std::to_string(t) + ". Try decreasing the learning rate");

        } else if(p_dat.U.has_inf() || p_dat.U.has_nan()){
          Rcpp::stop("ddhazard_fit_cpp estimation error: information matrix in correction step had inf or nan elements in bin " +
            std::to_string(t) + ". Try decreasing the learning rate");

        }

        if(p_dat.debug){
          Rcpp::Rcout << "Score vector and diagonal of information matrix at time " << t << " are:"<< std::endl;
          my_print(p_dat.u, "u");
          my_print(p_dat.U.diag(), "U");
        }

#ifdef USE_OPEN_BLAS
        openblas_set_num_threads(p_dat.n_threads);
#endif

        // E-step: scoring step: update values
        static bool inv_V_t_less_s_inv_have_failed;
        if(!arma::inv_sympd(V_t_less_s_inv, p_dat.V_t_less_s.slice(t - 1))){
          if(!inv_V_t_less_s_inv_have_failed){
            Rcpp::warning("V_(t|t-1) seems non-positive definite (or even worse non symmetric) at least once. Using general inverse instead");
            inv_V_t_less_s_inv_have_failed = true;
          }

          if(!arma::inv(V_t_less_s_inv, p_dat.V_t_less_s.slice(t - 1))){
            Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert V_(t|t-1)");
          }
        }

        static bool inv_V_t_t_s_inv_have_failed;
        arma::mat tmp_mat; // defined to avoid unhandled error if the next code throws
        if(!arma::inv_sympd(tmp_mat, V_t_less_s_inv + p_dat.U)){
          if(!inv_V_t_t_s_inv_have_failed){
            Rcpp::warning("V_(t|t) seems non-positive definite (or even worse non symmetric) at least once. Using general inverse instead");
            inv_V_t_t_s_inv_have_failed = true;
          }

          if(!arma::inv(tmp_mat, V_t_less_s_inv + p_dat.U)){
            Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to compute inverse for V_(t|t)");
          }
        }
        p_dat.V_t_t_s.slice(t) = std::move(tmp_mat); // arma objects are moveable

        p_dat.a_t_t_s.col(t) = i_a_t + p_dat.LR * p_dat.V_t_t_s.slice(t) * p_dat.u;

        if(!p_dat.is_mult_NR || arma::norm(p_dat.a_t_t_s.col(t) - i_a_t, 2) / (arma::norm(i_a_t, 2) + 1e-8) < p_dat.NR_eps)
          break;

        if(n_NR_it > p_dat.NR_it_max)
          Rcpp::stop("Failed to convergece in NR method of filter step within " + std::to_string(p_dat.NR_it_max) + " iterations");

        if(p_dat.debug){
          Rcpp::Rcout << "Did not converge in filter step in iteration " << n_NR_it << ". Convergence criteria value is  "
                      << arma::norm(p_dat.a_t_t_s.col(t) - i_a_t, 2) / (arma::norm(i_a_t, 2) + 1e-8) << std::endl;
        }

        i_a_t = p_dat.a_t_t_s.col(t);
      }

      p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * V_t_less_s_inv;

      if(p_dat.debug){
        std::stringstream str;
        str << t << "|" << t;

        my_print(p_dat.a_t_t_s.col(t), "a_(" + str.str() + ")");
        my_print(p_dat.V_t_t_s.slice(t).diag(), "diag(V_(" + str.str() + "))\n");
      }

      if(t == p_dat.d){
        arma::mat tmp_inv_mat;
        if(!arma::inv(tmp_inv_mat, arma::eye<arma::mat>(size(p_dat.U)) + p_dat.U * p_dat.V_t_less_s.slice(t - 1))){
          Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate for K_d matrix");
        }

        p_dat.K_d = p_dat.V_t_less_s.slice(t - 1) * tmp_inv_mat * p_dat.z_dot * diagmat(p_dat.H_diag_inv);
        // Parenthesis is key here to avoid making a n x n matrix for large n
        p_dat.K_d = (p_dat.F_ * p_dat.V_t_less_s.slice(t - 1) * p_dat.z_dot * diagmat(p_dat.H_diag_inv) * p_dat.z_dot.t()) * p_dat.K_d;
        p_dat.K_d = p_dat.F_ * p_dat.V_t_less_s.slice(t - 1) * p_dat.z_dot * diagmat(p_dat.H_diag_inv) -  p_dat.K_d;

        p_dat.lag_one_cov.slice(t - 1) = (arma::eye<arma::mat>(size(p_dat.U)) - p_dat.K_d * p_dat.z_dot.t()) * p_dat.F_ * p_dat.V_t_t_s.slice(t - 1);
      }
    }
  }
};





// This is the orginal UKF formulation from:
// Julier, Simon J., and Jeffrey K. Uhlmann. "New extension of the Kalman filter
// to nonlinear systems." AeroSense'97. International Society for Optics and
// Photonics, 1997

// Sigma points are though re-generated as suggested in:
// Menegaz, Henrique Marra Taira. "Unscented kalman filtering on euclidean and
// riemannian manifolds." (2016).

// The methods requires inversion of matrix with dimension equal to the dim of
// observational equation. Hence, it does not scale well in the number of
// observation per bin. The code is kept to test against
class UKF_solver_Org : public Solver{
  problem_data &p_dat;
  const uword m;
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
    for(uword i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2);
  }

public:
  UKF_solver_Org(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &k_):
  p_dat(p_),
  m(p_.a_t_t_s.n_rows),
  k(!k_.isNull() ? Rcpp::as< Rcpp::NumericVector >(k_)[0] : 3.0 - m),
  w_0(k / (m + k)),
  w_i(1 / (2 * (m + k))),
  sqrt_m_k(std::sqrt(m + k)),
  sigma_points(arma::mat(m, 2 * m + 1))
  {}

  void solve(){
#ifdef USE_OPEN_BLAS
    const int prev_n_thread = openblas_get_num_threads();
    openblas_set_num_threads(p_dat.n_threads);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif

    double event_time = p_dat.min_start;
    for (int t = 1; t < p_dat.d + 1; t++){

      double delta_t = p_dat.I_len[t - 1];
      event_time += delta_t;

      // Update sigma pooints
      compute_sigma_points(p_dat.a_t_t_s.unsafe_col(t - 1),
                           sigma_points, p_dat.V_t_t_s.slice(t - 1));

      // E-step: Filter step
      // First we compute the mean
      p_dat.a_t_less_s.col(t - 1) = w_0 * sigma_points.unsafe_col(0) +
        w_i * arma::sum(sigma_points.cols(1, sigma_points.n_cols - 1), 1);

      // Then the variance
      p_dat.V_t_less_s.slice(t - 1) = delta_t * p_dat.Q;
      for(uword i = 0; i < sigma_points.n_cols; ++i){
        const double &w = i == 0 ? w_0 : w_i;

        p_dat.V_t_less_s.slice(t - 1) +=
          (w * (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
          (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1)).t();
      }

      // Re-generate
      compute_sigma_points(p_dat.a_t_less_s.col(t - 1),
                           sigma_points, p_dat.V_t_less_s.slice(t - 1));

      // E-step: correction-step
      // First, we compute the mean of the outcome
      arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;

      arma::mat Z_t = (sigma_points.t() * p_dat.X.cols(r_set)).t(); // we transpose due to the column-major
      Z_t = arma::trunc_exp(Z_t);
      Z_t.for_each([](arma::mat::elem_type &val) { val = val / (1 + val); });

      // Compute y_bar, P_a_v and P_v_v
      const arma::vec y_bar = w_0 * Z_t.unsafe_col(0) +
        w_i * arma::sum(Z_t.cols(1, Z_t.n_cols - 1), 1);

      arma::mat P_a_v = (w_0 * (sigma_points.unsafe_col(0) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
        (Z_t.unsafe_col(0) - y_bar).t();

      arma::mat P_v_v = (w_0 * (Z_t.unsafe_col(0) - y_bar)) * (Z_t.unsafe_col(0) - y_bar).t() +
        arma::diagmat(w_0 * Z_t.unsafe_col(0) % (1 - Z_t.unsafe_col(0)));

      for(uword i = 1; i < sigma_points.n_cols; ++i){
        P_a_v += (w_i * (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
          (Z_t.unsafe_col(i) - y_bar).t();
        P_v_v += (w_i * (Z_t.unsafe_col(i) - y_bar)) * (Z_t.unsafe_col(i) - y_bar).t() +
          arma::diagmat(w_i * Z_t.unsafe_col(i) % (1 - Z_t.unsafe_col(i)));
      }

      P_v_v.diag() += p_dat.ridge_eps;

      // Compute new estimates
      if(!arma::inv_sympd(P_v_v, P_v_v)){ // NB: Note that we invert the matrix here so P_v_v is inv(P_v_v)
        Rcpp::warning("Failed to use inversion for symmetric square matrix for P_v_v. Trying with general inversion method");
        if(!arma::inv_sympd(P_v_v, P_v_v)){
          Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert P_v_v");
        }
      }

      p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.unsafe_col(t - 1) +
        P_a_v * (P_v_v * ((p_dat.is_event_in_bin(r_set) == t - 1) - y_bar));

      p_dat.V_t_t_s.slice(t) = p_dat.V_t_less_s.slice(t - 1) - P_a_v * P_v_v * P_a_v.t();

      p_dat.B_s.slice(t - 1) = arma::solve(p_dat.V_t_less_s.slice(t - 1), p_dat.F_ * p_dat.V_t_t_s.slice(t - 1)).t();
    }

#ifdef USE_OPEN_BLAS
    openblas_set_num_threads(prev_n_thread);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif
  }
};


// New method that use the Woodbury matrix identity to make the UKF algorithm
// scale linearly with the dimension of the observationally.

// beta = 0.0 and alpha = 1.0 yields the same sigma points as in:
// Julier, Simon J., and Jeffrey K. Uhlmann. "New extension of the Kalman filter
// to nonlinear systems." AeroSense'97. International Society for Optics and
// Photonics, 1997

// Altering these will yields parameter estimate similar to:
// Wan, Eric A., and Rudolph Van Der Merwe. "The unscented Kalman filter for
// nonlinear estimation." Adaptive Systems for Signal Processing, Communications
// and Control Symposium 2000. AS-SPCC. The IEEE 2000. Ieee, 2000.

// We re-generate the sigma points as suggested in:
// Menegaz, Henrique Marra Taira. "Unscented kalman filtering on euclidean and
// riemannian manifolds." (2016).
class UKF_solver_New : public Solver{
protected:
  problem_data &p_dat;
  const uword m;
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

  void compute_sigma_points(const arma::vec &a_t,
                            arma::mat &s_points,
                            const arma::mat &P_x_x){
    arma::mat cholesky_decomp;
    if(!arma::chol(cholesky_decomp, P_x_x, "lower")){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Cholesky decomposition failed");
    }

    s_points.col(0) = a_t;
    for(uword i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_lambda * cholesky_decomp.unsafe_col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_lambda * cholesky_decomp.unsafe_col((i - 1) / 2);
  }

public:
  UKF_solver_New(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                 Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                 Rcpp::Nullable<Rcpp::NumericVector> &beta):
  p_dat(p_),
  m(p_.a_t_t_s.n_rows),

  a(!alpha.isNull() ? Rcpp::as< Rcpp::NumericVector >(alpha)[0] : 1e-1),
  k(!kappa.isNull() ? Rcpp::as< Rcpp::NumericVector >(kappa)[0] : m / pow(a, 2.0) - m),
  b(!beta.isNull() ? Rcpp::as< Rcpp::NumericVector >(beta)[0] : 2.0),
  lambda(pow(a, 2) * (m + k) - m),

  w_0(lambda / (m + lambda)),
  w_0_c(w_0 + 1 - pow(a, 2) + b),
  w_0_cc(w_0 + 1 - a),
  w_i(1 / (2 * (m + lambda))),
  sqrt_m_lambda(std::sqrt(m + lambda)),

  sigma_points(arma::mat(m, 2 * m + 1))
  {
    weights_vec = arma::vec(std::vector<double>(2 * m + 1, w_i));
    weights_vec[0] = w_0;
    weights_vec_inv = arma::pow(weights_vec, -1);

    weights_vec_c = weights_vec;
    weights_vec_c[0] = w_0_c;
    weights_vec_c_inv = arma::pow(weights_vec_c, -1);

    weights_vec_cc = weights_vec;
    weights_vec_cc[0] = w_0_cc;

    if(p_dat.debug){
      Rcpp::Rcout << "alpha, beta, kappa, n_dim = "
                  << a << ", "
                  << b << ", "
                  << k << ", "
                  << m << std::endl;
      Rcpp::Rcout << "w_0, w_0_c, w_0_cc, w_i, lambda = "
                  << w_0 << ", "
                  << w_0_c << ", "
                  << w_0_cc << ", "
                  << w_i << ", "
                  << lambda << std::endl;
    }
  }

  void solve(){
#ifdef USE_OPEN_BLAS
    const int prev_n_thread = openblas_get_num_threads();
    openblas_set_num_threads(p_dat.n_threads);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif

    const arma::vec offsets = p_dat.any_fixed_in_M_step ?
      p_dat.fixed_terms.t() * p_dat.fixed_parems : arma::vec(p_dat.X.n_cols, arma::fill::zeros);
    double bin_stop = p_dat.min_start;
    for (int t = 1; t < p_dat.d + 1; t++){
      double bin_start = bin_stop;
      double delta_t = p_dat.I_len[t - 1];
      bin_stop += delta_t;

      // E-step: Filter step
      p_dat.a_t_less_s.col(t - 1) = p_dat.F_ *  p_dat.a_t_t_s.unsafe_col(t - 1);
      p_dat.V_t_less_s.slice(t - 1) = p_dat.F_ * p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ + delta_t * p_dat.Q;

      // Re-generate
      compute_sigma_points(p_dat.a_t_less_s.col(t - 1),
                           sigma_points, p_dat.V_t_less_s.slice(t - 1));

      if(p_dat.debug){
        my_print(p_dat.V_t_less_s.slice(t - 1), "Chol decomposing for regenerations:");
        my_print(sigma_points, "new sigma points");
      }

      // E-step: correction-step
      arma::mat O;
      arma::vec c_vec;
      arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;

      Compute_intermediates(r_set, offsets(r_set), t, bin_start, bin_stop, c_vec, O);

      // Substract mean to get delta sigma points
      arma::mat delta_sigma_points = sigma_points.each_col() - p_dat.a_t_less_s.unsafe_col(t - 1);

      p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.unsafe_col(t - 1) + p_dat.LR * delta_sigma_points * (weights_vec_cc % c_vec);

      p_dat.V_t_t_s.slice(t) = p_dat.V_t_less_s.slice(t - 1) -
        (delta_sigma_points.each_row() % weights_vec_cc.t()) * O * (delta_sigma_points.each_row() % weights_vec_cc.t()).t();

      if(p_dat.debug){
        std::stringstream str, str_less;
        str_less << t << "|" << t - 1;
        str << t << "|" << t;

        my_print(p_dat.V_t_less_s.slice(t - 1).diag(), "diag(V_(" + str_less.str() + "))");
        my_print(p_dat.V_t_t_s.slice(t).diag(), "diag(V_(" + str.str()  + "))");
        my_print(p_dat.a_t_less_s.col(t - 1), "a_(" + str_less.str() + ")");
        my_print(p_dat.a_t_t_s.col(t), "a_(" + str.str() + ")");
      }

      // Solve should be faster than inv_sympd http://arma.sourceforge.net/docs.html#inv_sympd
      // Solves yields a solution X to A * X = B <=> X = A^-1 B
      // We are looking at:
      //  X = B A^-1
      // X^T = A^-1 B^T <=> A X^T = B^T
      p_dat.B_s.slice(t - 1) = arma::solve(
        p_dat.V_t_less_s.slice(t - 1), p_dat.F_ * p_dat.V_t_t_s.slice(t - 1)).t();
    }

#ifdef USE_OPEN_BLAS
    openblas_set_num_threads(prev_n_thread);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif
  }
};


class UKF_solver_New_logit : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1: Compute expected outcomes given sigma points **
    O = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // we transpose due to the column-major

    O.each_col() += offsets;

    O = arma::trunc_exp(O);

    O.for_each([](arma::mat::elem_type &val) { val = val / (1 + val); });

    // ** 2: Compute mean observation sing sigma points **
    const arma::vec y_bar = w_0 * O.unsafe_col(0) +
      w_i * arma::sum(O.cols(1, O.n_cols - 1), 1);

    // ** 3: Compute variances and expected variance **
    arma::vec vars = w_0_c * (O.unsafe_col(0) % (1.0 - O.unsafe_col(0)));
    for(uword i = 1; i < O.n_cols; ++i){
      vars += w_i * (O.unsafe_col(i) % (1.0 - O.unsafe_col(i)));
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() / vars).t() * ((p_dat.is_event_in_bin(r_set) == t - 1) - y_bar);
    O = O.t() * (O.each_col() / vars);

    // Compute intermediate matrix
    arma::mat tmp_mat;
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_inv) + O)){ // this is symetric but not gauranteed to be postive definie due to ponetial negative weigths in weights_vec_inv
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O)){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_logit(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                       Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                       Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};


class UKF_solver_New_exp_bin : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk, sigma_points.n_cols);
    arma::vec vars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk length
    arma::vec at_risk_length(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      at_risk_length(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(-1 * v);

        O(j, i) = exp_model_funcs::expect_chance_die(v, inv_exp_v);
        vars(j) += w_c * exp_model_funcs::var_chance_die(v, inv_exp_v);
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() / vars).t() * (do_die - y_bar);
    O = O.t() * (O.each_col() / vars);

    // Compute intermediate matrix
    arma::mat tmp_mat;
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_inv) + O)){ // this is symmetric but not gauranteed to be postive definie due to ponetial negative weigths in weights_vec_inv
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O)){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exp_bin(
    problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
    Rcpp::Nullable<Rcpp::NumericVector> &alpha,
    Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};

class UKF_solver_New_exp_trunc_time : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk, sigma_points.n_cols);
    arma::vec vars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk length
    arma::vec at_risk_length(n_risk), time_outcome(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      time_outcome(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
      at_risk_length(i) = do_die(i) ?
      bin_tstop - std::max(p_dat.tstart(*it), bin_tstart) : time_outcome(i);
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(-1 * v);

        O(j, i) = exp_model_funcs::expect_time(v, at_risk_length(j), inv_exp_v, exp_eta);
        vars(j) += w_c * exp_model_funcs::var_wait_time(v, at_risk_length(j),  exp_eta, inv_exp_v);
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() / vars).t() * (time_outcome - y_bar);
    O = O.t() * (O.each_col() / vars);

    // Compute intermediate matrix
    arma::mat tmp_mat;
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_inv) + O)){ // this is symmetric but not gauranteed to be postive definie due to ponetial negative weigths in weights_vec_inv
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O)){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exp_trunc_time(
    problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
    Rcpp::Nullable<Rcpp::NumericVector> &alpha,
    Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};


class UKF_solver_New_exp_trunc_time_w_jump : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O){
    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk, sigma_points.n_cols);
    arma::vec vars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk length
    arma::vec at_risk_length(n_risk), time_outcome(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      time_outcome(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
      at_risk_length(i) = do_die(i) ?
        bin_tstop - std::max(p_dat.tstart(*it), bin_tstart) : time_outcome(i);

      if(do_die(i)){
        // we deduct the at risk lenght if the indvidual dies
        time_outcome(i) -= at_risk_length(i);
      }
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double inv_exp_eta = pow(exp_eta, -1);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(- v);

        O(j, i) = exp_model_funcs::expect_time_w_jump(exp_eta, inv_exp_eta, inv_exp_v, at_risk_length(j));
        vars(j) += w_c * exp_model_funcs::var_wait_time_w_jump(exp_eta, inv_exp_v, at_risk_length(j));
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Substract y_bar to get deviations
    O.each_col() -= y_bar;

    c_vec = (O.each_col() / vars).t() * (time_outcome - y_bar);
    O = O.t() * (O.each_col() / vars);

    // Compute intermediate matrix
    arma::mat tmp_mat;
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_inv) + O)){ // this is symmetric but not gauranteed to be postive definie due to ponetial negative weigths in weights_vec_inv
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // Compute vector for state space vector
    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    // Re-compute intermediate matrix using the other weight vector
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O)){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    // compute matrix for co-variance
    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exp_trunc_time_w_jump(
    problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
    Rcpp::Nullable<Rcpp::NumericVector> &alpha,
    Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {}
};


class UKF_solver_New_exponential : public UKF_solver_New{
  void Compute_intermediates(const arma::uvec &r_set,
                             const arma::vec offsets,
                             const int t,
                             const double bin_tstart, const double bin_tstop,
                             arma::vec &c_vec, arma::mat &O)
  {
    // See comments in UKF_solver_New_logit. The main difference here is that
    // we have tuples as outcomes. Thus, we have to deal with covariance terms

    // ** 1-3: compute outcome given sigma points, means and variances **
    const arma::uword n_risk = r_set.n_elem;
    O.set_size(n_risk * 2, sigma_points.n_cols);
    arma::vec vars(n_risk * 2, arma::fill::zeros);
    arma::vec covars(n_risk, arma::fill::zeros);

    arma::vec y_bar(n_risk * 2, arma::fill::zeros);

    arma::mat etas = (sigma_points.rows(p_dat.span_current_cov).t() * p_dat.X.cols(r_set)).t(); // linear predictors

    // Armadillo do not have a bool vector so we use an integer vector instead
    arma::ivec do_die = arma::conv_to<arma::ivec>::from(p_dat.is_event_in_bin(r_set) == t - 1);

    // We need two times: the length at risk and the outcome in case the
    // individual dies before the end of the bin
    arma::vec at_risk_length(n_risk), time_outcome(n_risk);

    auto it = r_set.begin();
    for(arma::uword i = 0; i < n_risk; ++i, ++it){
      time_outcome(i) = std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
      at_risk_length(i) = do_die(i) ?
      bin_tstop - std::max(p_dat.tstart(*it), bin_tstart) : time_outcome(i);
    }

    // Compute variance and mean
    for(arma::uword i = 0; i < sigma_points.n_cols; ++i){
      double w = (i == 0) ? w_0 : w_i;
      double w_c = (i == 0) ? w_0_c : w_i;

      const arma::vec exp_etas = arma::trunc_exp(offsets + etas.col(i));

      for(arma::uword j = 0; j < n_risk; ++j){
        const double exp_eta = exp_etas(j);
        const double v = at_risk_length(j) * exp_eta;

        const double inv_exp_v = exp(-1 * v);

        O(j, i) = exp_model_funcs::expect_chance_die(v, inv_exp_v);
        vars(j) += w_c * exp_model_funcs::var_chance_die(v, inv_exp_v);

        covars(j) += w_c * exp_model_funcs::covar(v, inv_exp_v, exp_eta);

        O(j + n_risk, i) = exp_model_funcs::expect_time(
          v, at_risk_length(j), inv_exp_v, exp_eta);
        vars(j + n_risk) += w_c * exp_model_funcs::var_wait_time(
          v, at_risk_length(j), exp_eta, inv_exp_v);
      }

      y_bar += w * O.col(i);
    }

    vars += p_dat.ridge_eps;

    // ** 4: Compute c **
    // Defines span to avoid 40-error
    arma::span span_binary(0, n_risk - 1);
    arma::span span_time(n_risk, 2 * n_risk -1);

    // Compute diagonal terms of inverse covariance matrix
    arma::vec inv_covmat_diag(2 * n_risk);
    inv_covmat_diag.subvec(span_binary) =
      (vars.subvec(span_binary) - covars % arma::pow(vars.subvec(span_time), -1) % covars);

    inv_covmat_diag.subvec(span_time) =
      (vars.subvec(span_time) - covars % arma::pow(vars.subvec(span_binary), -1) % covars);

    inv_covmat_diag.transform([](double val) { return(1 / val); });

    // Compute off diagonal terms of inverse covariance matrix
    arma::vec inv_covmat_off_diag(n_risk);
    inv_covmat_off_diag = -1 * inv_covmat_diag.subvec(span_time) % covars / vars.subvec(span_binary);

    // Compute inverse covariance matrix dot centered sigma points
    O.each_col() -= y_bar; // center

    O = O.t(); // transpose to ease the next computations

    arma::mat tmp_mat = O.each_row() % inv_covmat_diag.t();
    tmp_mat.cols(span_binary) +=
      O.cols(span_time).each_row() %  inv_covmat_off_diag.t();
    tmp_mat.cols(span_time) +=
      O.cols(span_binary).each_row() %  inv_covmat_off_diag.t();

    O = O.t(); // transpose back

    {
      arma::vec outcome(n_risk * 2);
      outcome.subvec(span_binary) = arma::conv_to<arma::vec>::from(do_die);
      outcome.subvec(span_time) = time_outcome;

      c_vec = tmp_mat * (outcome - y_bar);
    }

    O = tmp_mat * O;
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_inv) + O)){ // this is symetric but not gauranteed to be postive definie due to ponetial negative weigths in weights_vec_inv
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    c_vec = c_vec -  tmp_mat * c_vec;

    // ** 5: Compute L using the notation in vignette **
    if(!arma::inv(tmp_mat, arma::diagmat(weights_vec_c_inv) + O)){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Failed to invert intermediate matrix in the scoring step");
    }
    tmp_mat = O * tmp_mat;

    O = O - tmp_mat * O;
  }

public:
  UKF_solver_New_exponential(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &kappa,
                             Rcpp::Nullable<Rcpp::NumericVector> &alpha,
                             Rcpp::Nullable<Rcpp::NumericVector> &beta):
  UKF_solver_New(p_, kappa, alpha, beta)
  {
  }
};















extern std::vector<double> logLike_cpp(const arma::mat&, const Rcpp::List&,
                                       const arma::mat&, const arma::mat&,
                                       arma::mat Q, const arma::mat&,
                                       const arma::vec&, const arma::vec&,
                                       const arma::vec &,
                                       const int, const std::string);



// Method to estimate fixed effects like in biglm::bigglm
template<typename T>
void estimate_fixed_effects(problem_data * const p_data, const int chunk_size,
                            bigglm_updateQR<T> &updater){
  int it_outer = 0;
  arma::vec old_beta;
  do{
    int cursor_risk_set = 0;
    int n_elements = 0;

    // Set up look variables
    int t = 1; // start looping at one to be consistent with other implementations
    arma::mat fixed_terms(p_data->fixed_parems.n_elem, chunk_size, arma::fill::none);
    arma::vec offsets(chunk_size, arma::fill::zeros);
    arma::vec y(chunk_size);
    arma::vec eta;
    qr_obj qr(p_data->fixed_parems.n_elem);
    auto it = p_data->risk_sets.begin();
    double bin_stop = p_data->min_start;

    for(; it != p_data->risk_sets.end(); ++it, ++t){

      double bin_start = bin_stop;
      double delta_t = p_data->I_len[t - 1];
      bin_stop += delta_t;

      // Find the risk set and the number of elements to take
      arma::uvec r_set = Rcpp::as<arma::uvec>(*it) - 1;
      int r_set_size = r_set.n_elem;
      int n_elements_to_take = std::min(chunk_size - n_elements, r_set_size - cursor_risk_set);
      r_set = r_set.subvec(cursor_risk_set, cursor_risk_set + n_elements_to_take - 1);

      // Find the outcomes, fixed terms and compute the offsets
      y.subvec(n_elements, n_elements + n_elements_to_take - 1) =
        arma::conv_to<arma::vec>::from(p_data->is_event_in_bin.elem(r_set) == (t - 1));

      fixed_terms.cols(n_elements, n_elements + n_elements_to_take - 1) =
        p_data->fixed_terms.cols(r_set);

      if(p_data->any_dynamic){
        offsets.subvec(n_elements, n_elements + n_elements_to_take - 1) =
          p_data->X.cols(r_set).t() * p_data->a_t_t_s.col(t).head(p_data->n_params_state_vec);
      } else {
        offsets.subvec(n_elements, n_elements + n_elements_to_take - 1).fill(0.);
      }

      for(arma::uword i = 0; i < r_set.n_elem; ++i){
        offsets(n_elements + i) +=
          T().time_offset(std::min(p_data->tstop(r_set(i)), bin_stop)
                            - std::max(p_data->tstart(r_set(i)), bin_start));
      }

      n_elements += n_elements_to_take;

      if(n_elements == chunk_size){ // we have reached the chunk_size

        arma::vec eta = fixed_terms.t() * p_data->fixed_parems;
        updater.update(qr, fixed_terms, eta, offsets, y);

        n_elements = 0;
      } else if(it == --p_data->risk_sets.end()){ // there is no more bins to process

        y = y.subvec(0, n_elements - 1);
        fixed_terms = fixed_terms.cols(0, n_elements - 1);
        offsets = offsets.subvec(0, n_elements - 1);

        arma::vec eta =  fixed_terms.t() * p_data->fixed_parems;
        updater.update(qr, fixed_terms, eta, offsets, y);
      }

      if(cursor_risk_set + n_elements_to_take < r_set_size){ // there are still elements left in the bin
        cursor_risk_set = cursor_risk_set + n_elements_to_take;
        --it;
        --t;
        bin_stop -= delta_t;
      } else
        cursor_risk_set = 0;
    }

    old_beta = p_data->fixed_parems;
    p_data->fixed_parems = bigglm_regcf(qr);

  } while(++it_outer < p_data->max_it_fixed_params && // Key that this the first condition we check when we use &&
    arma::norm(p_data->fixed_parems - old_beta, 2) / (arma::norm(old_beta, 2) + 1e-8) > p_data->eps_fixed_parems);

  static bool failed_to_converge_once = false;
  if(it_outer == p_data->max_it_fixed_params && !failed_to_converge_once){
    failed_to_converge_once = true;
    std::stringstream msg;
    msg << "Failed to estimate fixed effects in " << p_data->max_it_fixed_params << " iterations at least once" << std::endl;
    Rcpp::warning(msg.str());
  }
}



// [[Rcpp::export]]
Rcpp::List ddhazard_fit_cpp(arma::mat &X, arma::mat &fixed_terms, // Key: assumed to have observations in the columns for performance due to column-major storage
                            const arma::vec &tstart, const arma::vec &tstop,
                            const arma::colvec &a_0,
                            const arma::vec &fixed_parems_start,
                            arma::mat Q_0, // by value copy. This  is key cuz we will change it if est_Q_0 = T
                            arma::mat Q, // similary this is a copy
                            const Rcpp::List &risk_obj,
                            const arma::mat &F_,
                            const double eps_fixed_parems, const int max_it_fixed_params,
                            const arma::uword n_max = 100, const double eps = 0.001,
                            const arma::uword verbose = 0,
                            const int order_ = 1, const bool est_Q_0 = true,
                            const std::string method = "EKF",
                            Rcpp::Nullable<Rcpp::NumericVector> kappa = R_NilValue, // see this link for nullable example http://blogs.candoerz.com/question/164706/rcpp-function-for-adding-elements-of-a-vector.aspx
                            Rcpp::Nullable<Rcpp::NumericVector> alpha = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> beta = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> NR_eps = R_NilValue,
                            Rcpp::Nullable<Rcpp::NumericVector> LR = R_NilValue,
                            const std::string model = "logit",
                            const std::string M_step_formulation = "Fahrmier94",
                            const int fixed_effect_chunk_size = 2e4,
                            const bool debug = false,
                            const unsigned int NR_it_max = 100,
                            const int n_threads = -1,
                            const double ridge_eps = .0001,
                            const int n_fixed_terms_in_state_vec = 0){
  if(Rcpp::as<bool>(risk_obj["is_for_discrete_model"]) &&
     is_exponential_model(model)){
    Rcpp::stop("risk_obj has 'is_for_discrete_model' = true which should be false for model '" + model  +"'");
  } else if(!Rcpp::as<bool>(risk_obj["is_for_discrete_model"]) && model == "logit"){
    Rcpp::stop("risk_obj has 'is_for_discrete_model' = false which should be true for model '" + model  +"'");
  }

  // Declare non constants and intialize some of them
  double delta_t, test_max_diff;
  const double Q_warn_eps = sqrt(std::numeric_limits<double>::epsilon());

  Rcpp::NumericVector conv_values;

  uword it = 0;

  // M-stp pointers for convenience
  arma::mat *B, *V_less, *V;
  arma::vec a_less, a;

  const arma::ivec is_event_in_bin = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  // Intialize the solver for the E-step
  problem_data *p_data;
  Solver  *solver = NULL;

  if(method == "EKF"){
    p_data = new problem_data_EKF(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      NR_eps, LR,
      eps_fixed_parems, max_it_fixed_params,
      n_max, eps, verbose,
      order_, est_Q_0, model != "logit", NR_it_max, debug, n_threads,
      ridge_eps);
    solver = new EKF_solver(static_cast<problem_data_EKF &>(*p_data), model);

  } else if (method == "UKF"){
    if(model != "logit" &&
       !is_exponential_model(model))
      Rcpp::stop("UKF is not implemented for model '" + model  +"'");
    p_data = new problem_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      eps_fixed_parems, max_it_fixed_params,
      n_max, eps, verbose,
      order_, est_Q_0, debug, LR, n_threads, ridge_eps);

    if(model == "logit"){
      solver = new UKF_solver_New_logit(*p_data, kappa, alpha, beta);

    } else if (model == "exp_combined"){
      solver = new UKF_solver_New_exponential(*p_data, kappa, alpha, beta);

    } else if (model == "exp_bin"){
      solver = new UKF_solver_New_exp_bin(*p_data, kappa, alpha, beta);
    } else if (model == "exp_trunc_time"){
      solver = new UKF_solver_New_exp_trunc_time(*p_data, kappa, alpha, beta);
    } else if (model == "exp_trunc_time_w_jump"){
      solver = new UKF_solver_New_exp_trunc_time_w_jump(*p_data, kappa, alpha, beta);
    } else
      Rcpp::stop("Model '", model ,"' is not implemented with UKF");

  } else if (method == "UKF_org"){
    if(model != "logit")
      Rcpp::stop("UKF is not implemented for model '" + model  +"'");

    p_data = new problem_data(
      n_fixed_terms_in_state_vec,
      X, fixed_terms, tstart, tstop, is_event_in_bin,
      a_0, fixed_parems_start, Q_0, Q,
      risk_obj, F_,
      eps_fixed_parems, max_it_fixed_params,
      n_max, eps, verbose,
      order_, est_Q_0, debug);

    if(p_data->any_fixed_in_M_step)
      Rcpp::stop("Fixed effects is not implemented with UKF");

    solver = new UKF_solver_Org(*p_data, kappa);
  }else{
    Rcpp::stop("method '" + method  + "'is not implemented");
  }

  arma::mat a_prev;
  a_prev.copy_size(p_data->a_t_t_s);
  a_prev.zeros();

  do
  {
    if(p_data->debug){
      if(it > 0)
        Rcpp::Rcout << "\n\n\n";
      Rcpp::Rcout << "##########################################\nStarting iteration " << it
                  << " with the following values" << std::endl;
      my_print(p_data->a_t_t_s.col(0), "a_0");
      my_print(p_data->Q.diag(), "diag(Q)");
    }


    if((it + 1) % 25 == 0)
      Rcpp::checkUserInterrupt(); // this is expensive (on Windows) - you do not want to check too often


    if(p_data->any_dynamic){
      p_data->V_t_t_s.slice(0) = Q_0; // Q_0 may have been updated or not

      // E-step
      solver->solve();

      // E-step: smoothing
      if(p_data->debug){
        Rcpp::Rcout << "Started smoothing" << std::endl;
      }

      for (int t = p_data->d - 1; t > -1; t--){
        // we need to compute the correlation matrix first
        if(t > 0){
          p_data->lag_one_cov.slice(t - 1) = p_data->V_t_t_s.slice(t) * p_data->B_s.slice(t - 1).t() +
            p_data->B_s.slice(t) * (
                p_data->lag_one_cov.slice(t) - F_ * p_data->V_t_t_s.slice(t)) * p_data->B_s.slice(t - 1).t();
        }

        p_data->a_t_t_s.col(t) = p_data->a_t_t_s.unsafe_col(t) + p_data->B_s.slice(t) *
          (p_data->a_t_t_s.unsafe_col(t + 1) - p_data->a_t_less_s.unsafe_col(t));
        p_data->V_t_t_s.slice(t) = p_data->V_t_t_s.slice(t) + p_data->B_s.slice(t) *
          (p_data->V_t_t_s.slice(t + 1) - p_data->V_t_less_s.slice(t)) * p_data->B_s.slice(t).t();

        if(p_data->debug){
          std::stringstream ss;
          ss << t << "|" <<  p_data->d;
          my_print(p_data->a_t_t_s.col(t), "a(" + ss.str() + ")");
          my_print(p_data->V_t_t_s.slice(t).diag(), "diag(V(" + ss.str() + "))");
        }
      }

      // M-step
      if(est_Q_0){
        Q_0 = p_data->V_t_t_s.slice(0);
      }
      Q.zeros();
      for (int t = 1; t < p_data->d + 1; t++){
        delta_t = p_data->I_len[t - 1];

        V_less = &p_data->V_t_t_s.slice(t - 1);
        V = &p_data->V_t_t_s.slice(t);
        a_less = p_data->a_t_t_s.unsafe_col(t - 1);
        a = p_data->a_t_t_s.unsafe_col(t);

        if(M_step_formulation == "Fahrmier94"){
          B = &p_data->B_s.slice(t - 1);

          Q += ((a - F_ * a_less) * (a - F_ * a_less).t() + *V
                  - F_ * *B * *V
                  - (F_ * *B * *V).t()
                  + F_ * *V_less * p_data->T_F_) / delta_t;

        } else if (M_step_formulation == "SmoothedCov"){
          B = &p_data->lag_one_cov.slice(t - 1); // this is not B but the lagged one smooth correlation. Do not mind the variable name

          Q += ((a - F_ * a_less) * (a - F_ * a_less).t() + *V
                  - F_ * *B
                  - (F_ * *B).t()
                  + F_ * *V_less * p_data->T_F_) / delta_t;
        } else
          Rcpp::stop("'M_step_formulation' of type '" + M_step_formulation + "' is not implemented");

      }
      Q /= p_data->d;

      if(p_data->debug){
        my_print(p_data->Q.diag(), "Diag(Q) before changes in M-step");
      }


      if(p_data->any_fixed_terms_in_state_vec){
        Q.rows(p_data->span_fixed_params).zeros();
        Q.cols(p_data->span_fixed_params).zeros();
      }

      if((test_max_diff = static_cast<arma::mat>(Q - Q.t()).max()) > Q_warn_eps){
        std::ostringstream warning;
        warning << "Q - Q.t() maximal element difference was " << test_max_diff <<
          " in iteration " << it + 1;
        Rcpp::warning(warning.str());
      }

      if((test_max_diff = static_cast<arma::mat>(Q_0 - Q_0.t()).max()) > Q_warn_eps){
        std::ostringstream warning;
        warning << "Q_0 - Q_0.t() maximal element difference was " << test_max_diff <<
          " in iteration " << it + 1;
        Rcpp::warning(warning.str());
      }

      // Ensure that Q and Q_0 are symmetric
      Q = (Q + Q.t()) / 2.0;
      Q_0 = (Q_0 + Q_0.t()) / 2.0;

      if(order_ > 1){
        arma::mat tmp_Q = Q(p_data->span_current_cov, p_data->span_current_cov);
        Q.zeros();
        Q(p_data->span_current_cov , p_data->span_current_cov) = tmp_Q;
      }
    }

    if(p_data->debug){
      my_print(p_data->Q.diag(), "Diag(Q) after changes in M-step");
    }

    if(p_data->any_fixed_terms_in_state_vec ||
       p_data->any_dynamic){
      conv_values.push_back(conv_criteria(a_prev(p_data->span_current_cov, arma::span::all),
                                          p_data->a_t_t_s(p_data->span_current_cov, arma::span::all)));
    } else
      conv_values.push_back(0.0);

    if(p_data->any_fixed_in_M_step){
      arma::vec old = p_data->fixed_parems;

      if(model == "logit"){
        bigglm_updateQR_logit  updater;
        estimate_fixed_effects(p_data, fixed_effect_chunk_size, updater);

      } else if(is_exponential_model(model)){
        bigglm_updateQR_poisson updater;
        estimate_fixed_effects(p_data, fixed_effect_chunk_size, updater);

      } else
        Rcpp::stop("Fixed effects is not implemented for '" + model  +"'");

      *(conv_values.end() -1) += conv_criteria(old, p_data->fixed_parems);
    }

    if(!p_data->any_dynamic) // No reason to take further iterations
      break;

    if(verbose && it % 5 < verbose){
      auto rcout_width = Rcpp::Rcout.width();

      arma::mat varying_only_F = p_data->F_; // take copy
      arma::mat varying_only_a = p_data->a_t_t_s; // take copy
      arma::vec fixed_effects_offsets;

      if(p_data->any_fixed_in_M_step){
        fixed_effects_offsets = p_data->fixed_terms.t() * p_data->fixed_parems;

      } else if(p_data->any_fixed_terms_in_state_vec){
        fixed_effects_offsets =
          p_data->X(p_data->span_fixed_params, arma::span::all).t() *
            p_data->a_t_t_s(p_data->span_fixed_params, arma::span::all).col(0);

        varying_only_a.shed_rows(p_data->span_fixed_params.a, p_data->span_fixed_params.b);
        varying_only_F.shed_rows(p_data->span_fixed_params.a, p_data->span_fixed_params.b);
        varying_only_F.shed_cols(p_data->span_fixed_params.a, p_data->span_fixed_params.b);


      } else{
        fixed_effects_offsets = arma::vec(p_data->X.n_cols, arma::fill::zeros);

      }

      double log_like =
        logLike_cpp(p_data->X(p_data->span_current_cov_varying, arma::span::all),
                    risk_obj,
                    varying_only_F,
                    Q_0(p_data->span_current_cov_varying, p_data->span_current_cov_varying),
                    Q(p_data->span_current_cov_varying, p_data->span_current_cov_varying),
                    varying_only_a,
                    p_data->tstart, p_data->tstop,
                    fixed_effects_offsets, order_, model)[0];
      Rcpp::Rcout << "Iteration " <<  std::setw(5)<< it + 1 <<
        " ended with conv criteria " << std::setw(15) << *(conv_values.end() -1) <<
          "\t" << "The log likelihood is " << log_like <<
            std::setw(rcout_width) << std::endl;
    }

    if(*(conv_values.end() -1) < eps)
      break;

    a_prev = p_data->a_t_t_s;
  }while(++it < n_max);

  if(it == n_max)
    throw std::runtime_error("EM algorithm did not converge within the n_max number of iterations");

  return(Rcpp::List::create(Rcpp::Named("V_t_d_s") = Rcpp::wrap(p_data->V_t_t_s),
                            Rcpp::Named("a_t_d_s") = Rcpp::wrap(p_data->a_t_t_s.t()),
                            Rcpp::Named("B_s") = Rcpp::wrap(p_data->B_s),
                            Rcpp::Named("lag_one_cov") = Rcpp::wrap(p_data->lag_one_cov),
                            Rcpp::Named("fixed_effects") = Rcpp::wrap(p_data->fixed_parems),

                            Rcpp::Named("n_iter") = it + 1,
                            Rcpp::Named("conv_values") = conv_values,
                            Rcpp::Named("Q") = Rcpp::wrap(Q),
                            Rcpp::Named("Q_0") = Rcpp::wrap(Q_0)));
}
