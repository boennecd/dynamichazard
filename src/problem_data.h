#ifndef DDHAZARD_DATA
#define DDHAZARD_DATA

#include "arma_n_rcpp.h"
#include <future>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Object used for map function. We keep a unique pointer to a new object if it
 * is created in the function to make sure the orginal object is
 not destructed. */
template <typename T_view, typename T_type>
class map_res {
  using ptr_T = std::unique_ptr<T_type>;

public:
  T_view subview;
  ptr_T org_ptr;

  map_res(T_view subview): subview(subview) {}
  map_res(T_view subview, ptr_T &ptr): subview(subview) {
    org_ptr = std::move(ptr);
  }
};

using map_res_col = map_res<arma::subview_col<double>, arma::vec>;
using map_res_mat = map_res<arma::subview<double>, arma::mat>;

class problem_data {
public:
  // Constants
  const bool any_dynamic;         // Any time-varying coefficients?
  const bool any_fixed_in_E_step;
  const bool any_fixed_in_M_step;

  const int d;
  const Rcpp::List risk_sets;

  const unsigned int space_dim; // dimension of state space vector
  const unsigned int covar_dim; // dimension of dynamic covariate vector

  // TODO: only relevant for random walk models -- move to derived class
  const int order;
  const int n_params_state_vec_fixed;
  const int n_params_state_vec_varying; // NB: not including order
  const int n_params_state_vec;         // NB: not including order

  // indicies for dot product for linear predictor. Includes both time-varying
  // effects and fixed effects estimated in the E-step
  const std::unique_ptr<const arma::span> span_current_cov;
  // indicies for varying parameters in dot product
  const std::unique_ptr<const arma::span> span_current_cov_varying;
  // indicies of fixed terms in state vector
  const std::unique_ptr<const arma::span> span_fixed_params;

  // TODO: end of "only relevant ..."

  const arma::mat &F_;
  const arma::mat T_F_;

  arma::mat X;
  arma::mat fixed_terms; // used if fixed terms are estimated in the M-step

  const std::vector<double> I_len;

  const int n_threads;

  const arma::vec &tstart;
  const arma::vec &tstop;
  const arma::ivec &is_event_in_bin;

  const double min_start;

  // non-constants
  arma::mat &Q;
  arma::mat &Q_0;

  problem_data(const int n_fixed_terms_in_state_vec,
               arma::mat &X,
               arma::mat &fixed_terms,
               const arma::vec &tstart,
               const arma::vec &tstop, const arma::ivec &is_event_in_bin,
               const arma::colvec &a_0,
               arma::mat &Q_0,
               arma::mat &Q,
               const Rcpp::List &risk_obj,
               const arma::mat &F_,
               const int n_max,
               const int order,
               const int n_threads):
    any_dynamic(X.n_elem > 0),
    any_fixed_in_E_step(n_fixed_terms_in_state_vec > 0),
    any_fixed_in_M_step(fixed_terms.n_elem > 0),

    d(Rcpp::as<int>(risk_obj["d"])),
    risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"])),

    space_dim(a_0.size()),
    covar_dim(X.n_rows),

    order(order),
    n_params_state_vec_fixed(n_fixed_terms_in_state_vec),
    n_params_state_vec_varying((a_0.size() - n_fixed_terms_in_state_vec) / order),
    n_params_state_vec(n_params_state_vec_fixed + n_params_state_vec_varying),


    // Seems like there is no empty span for arma
    // See: armadillo-code/include/armadillo_bits/span.hpp
    //      and https://stackoverflow.com/q/38155952/5861244
    // Thus, I decided to use const std::unique_ptr<const> and set it to
    // nullptr when the span should not be used
    span_current_cov(
      (any_dynamic || any_fixed_in_E_step) ?
        new arma::span(0, n_params_state_vec - 1) : nullptr),
    span_current_cov_varying(
      any_dynamic ?
        new arma::span(0, n_params_state_vec_varying - 1) : nullptr),
    span_fixed_params(
      any_fixed_in_E_step ?
        new arma::span(n_params_state_vec_varying, n_params_state_vec - 1) : nullptr),

    F_(F_),
    T_F_((any_dynamic || any_fixed_in_E_step) ? F_.t() : arma::mat()),

    X(X.begin(), X.n_rows, X.n_cols, false),
    fixed_terms(fixed_terms.begin(), fixed_terms.n_rows, fixed_terms.n_cols, false),
    I_len(Rcpp::as<std::vector<double> >(risk_obj["I_len"])),

    n_threads((n_threads > 0) ? n_threads : std::thread::hardware_concurrency()),

    tstart(tstart),
    tstop(tstop),
    is_event_in_bin(is_event_in_bin),
    min_start(Rcpp::as<double>(risk_obj["min_start"])),

    Q(Q),
    Q_0(Q_0)
  {
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
    omp_set_nested(0);
#endif
  }

  problem_data& operator=(const problem_data&) = delete;
  problem_data(const problem_data&) = delete;
  problem_data() = delete;
};

/* problem_data used in ddhazard_fit.cpp */
class ddhazard_data : public problem_data {
public:
  // constants
  const arma::vec &weights;

  const double event_eps; // something small
  const double eps_fixed_parems;
  const int max_it_fixed_params;

  const double denom_term;

  const bool debug;
  const double LR;

  const bool use_pinv;
  const std::string criteria;

  // Declare non constants. Some are intialize
  arma::vec fixed_parems;
  arma::mat a_t_t_s;
  arma::mat a_t_less_s;

  arma::cube V_t_t_s;
  arma::cube V_t_less_s;
  arma::cube B_s;

  arma::cube lag_one_cov;

  // Information for debugging
  std::string computation_stage;
  int em_iteration;

  ddhazard_data(const int n_fixed_terms_in_state_vec,
                arma::mat &X,
                arma::mat &fixed_terms,
                const arma::vec &tstart,
                const arma::vec &tstop, const arma::ivec &is_event_in_bin,
                const arma::colvec &a_0,
                const arma::vec &fixed_parems_start,
                arma::mat &Q_0,
                arma::mat &Q,
                const Rcpp::List &risk_obj,
                const arma::mat &F_,
                const double eps_fixed_parems,
                const int max_it_fixed_params,
                const arma::vec &weights,
                const int n_max, const double eps,
                const bool verbose,
                const int order, const bool est_Q_0,
                const bool debug,
                Rcpp::Nullable<Rcpp::NumericVector> LR,
                const int n_threads,
                const double denom_term,
                const bool use_pinv,
                const std::string criteria):
    problem_data(
      n_fixed_terms_in_state_vec,
      X,
      fixed_terms,
      tstart,
      tstop, is_event_in_bin,
      a_0,
      Q_0,
      Q,
      risk_obj,
      F_,
      n_max,
      order,
      n_threads),
    weights(weights),

    event_eps(d * std::numeric_limits<double>::epsilon()),
    eps_fixed_parems(eps_fixed_parems),
    max_it_fixed_params(max_it_fixed_params),
    denom_term(denom_term),

    debug(debug),
    LR(LR.isNotNull() ? Rcpp::as< Rcpp::NumericVector >(LR)[0] : 1.0),

    use_pinv(use_pinv),
    criteria(criteria),
    fixed_parems(fixed_parems_start)
  {
    if(debug)
      Rcpp::Rcout << "Using " << n_threads << " threads" << std::endl;

    if(any_dynamic || any_fixed_in_E_step){
      a_t_t_s = arma::mat(space_dim, d + 1);
      a_t_less_s = arma::mat(space_dim, d);
      V_t_t_s = arma::cube(space_dim, space_dim, d + 1);
      V_t_less_s = arma::cube(space_dim, space_dim, d);
      B_s = arma::cube(space_dim, space_dim, d);

      a_t_t_s.col(0) = a_0;

      lag_one_cov = arma::cube(space_dim, space_dim, d);
    }
  }

  ddhazard_data & operator=(const ddhazard_data&) = delete;
  ddhazard_data(const ddhazard_data&) = delete;
  ddhazard_data() = delete;

  map_res_col get_dynamic_coefs(unsigned int);

  /* maps from state space dimension to covariate dimension */
  virtual map_res_col lp_map(arma::vec&) = 0;
  virtual map_res_mat lp_map(arma::mat&) = 0;

  /* inverse of the above */
  virtual map_res_col lp_map_inv(arma::vec&) = 0;
  virtual map_res_mat lp_map_inv(arma::mat&) = 0;
};

inline std::string debug_msg_prefix(const ddhazard_data &dat){
  std::stringstream out;
  out << "--it " << std::setw(5) <<  dat.em_iteration
      << ", " << dat.computation_stage << ": ";
  return(out.str());
}

template<typename T>
inline
void
my_print(const ddhazard_data &dat, const T &X, std::string msg = "")
{
  my_print(X, msg, debug_msg_prefix(dat));
}

class my_debug_logger{
private:
  const ddhazard_data *dat;

protected:
  std::ostringstream os;

public:
  my_debug_logger(const ddhazard_data &dat_):
  dat(&dat_) {}

  template<typename T>
  std::ostringstream& operator<<(const T &obj);

  ~my_debug_logger(){
    os << std::endl;
    Rcpp::Rcout << os.str();
  }
};

template<typename T>
std::ostringstream& my_debug_logger::operator<<(const T &obj){
  os << debug_msg_prefix(*dat) << obj;
  return os;
}

/* Concrete implementation of data class */
class ddhazard_data_random_walk : public ddhazard_data {
public:
  using ddhazard_data::ddhazard_data;

  map_res_col lp_map(arma::vec&);
  map_res_mat lp_map(arma::mat&);
  map_res_col lp_map_inv(arma::vec&);
  map_res_mat lp_map_inv(arma::mat&);
};

#endif
