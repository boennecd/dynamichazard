#ifndef DDHAZARD_DATA
#define DDHAZARD_DATA

#include "arma_n_rcpp.h"
#include <future>
#include <type_traits>

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
  T_view sv; /* sv for (s)ub(v)iew */
  ptr_T org_ptr;

  map_res(T_view sv): sv(sv) {}
  map_res(T_view sv, ptr_T &ptr): sv(sv) {
    org_ptr = std::move(ptr);
  }
};

using map_res_col = map_res<arma::subview_col<double>, arma::vec>;
using map_res_mat = map_res<arma::subview<double>, arma::mat>;

enum side { left, both, right };

class problem_data {
protected:
  const arma::mat &F_;
  const arma::mat &R;
  const arma::mat &L;
  const arma::vec &m;

  using ptr_vec = std::unique_ptr<arma::vec>;
  using ptr_mat = std::unique_ptr<arma::mat>;

public:
  // Constants
  const bool any_dynamic;         // Any time-varying coefficients?
  const bool any_fixed_in_E_step;
  const bool any_fixed_in_M_step;

  const int d;
  const Rcpp::List risk_sets;

  const unsigned int space_dim; // dimension of state space vector
  const unsigned int covar_dim; // dimension of dynamic covariate vector
  const int n_params_state_vec_fixed; // TODO: maybe move to derived class

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

  problem_data(
    const int n_fixed_terms_in_state_vec,
    arma::mat &X,
    arma::mat &fixed_terms,
    const arma::vec &tstart,
    const arma::vec &tstop, const arma::ivec &is_event_in_bin,
    const arma::colvec &a_0,
    const arma::mat &R,
    const arma::mat &L,
    const arma::vec &m,
    arma::mat &Q_0,
    arma::mat &Q,
    const Rcpp::List &risk_obj,
    const arma::mat &F_,
    const int n_max,
    const int n_threads) :
    F_(F_), R(R), L(L), m(m),

    any_dynamic(X.n_elem > 0),
    any_fixed_in_E_step(n_fixed_terms_in_state_vec > 0),
    any_fixed_in_M_step(fixed_terms.n_elem > 0),

    d(Rcpp::as<int>(risk_obj["d"])),
    risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"])),

    space_dim(a_0.size()),
    covar_dim(X.n_rows),
    n_params_state_vec_fixed(n_fixed_terms_in_state_vec),

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

  /* maps previous to next state */
  map_res_col state_trans_map(arma::subview_col<double> a){
    arma::vec tmp(a.colptr(0), a.n_elem, false);
    return state_trans_map(tmp);
  }
  virtual map_res_col state_trans_map(arma::vec &a){
    ptr_vec ptr(new arma::vec(F_ * a + m));
    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

  virtual map_res_mat state_trans_map(arma::mat &M, side s = both){
    ptr_mat ptr;

    switch(s){
    case  left:
      ptr.reset(new arma::mat(F_ * M));
      break;
    case both:
      ptr.reset(new arma::mat(F_ * M * F_.t()));
      break;
    case right:
      ptr.reset(new arma::mat(     M * F_.t()));
      break;
    default:
      Rcpp::stop("'Side' not implemented");
    }

    arma::mat &out = *ptr.get();
    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }

  /* maps from errors to state space dimension */
  map_res_col err_state_map(arma::subview_col<double> a){
    arma::vec tmp(a.colptr(0), a.n_elem, false);
    return err_state_map(tmp);
  }
  virtual map_res_col err_state_map(arma::vec &a){
    ptr_vec ptr(new arma::vec(R * a));
    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

  virtual map_res_mat err_state_map(arma::mat &M){
    ptr_mat ptr(new arma::mat(R * M * R.t()));

    arma::mat &out = *ptr.get();
    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }

  /* inverse of the above */
  virtual map_res_mat err_state_map_inv(arma::mat &M){
   ptr_mat ptr(new arma::mat(R.t() * M * R));

     arma::mat &out = *ptr.get();
     return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }

  /* maps from state space dimension to covariate dimension */
  map_res_col lp_map(arma::subview_col<double> a){
    arma::vec tmp(a.colptr(0), a.n_elem, false);
    return lp_map(tmp);
  }
  virtual map_res_col lp_map(arma::vec&) = 0;
  virtual map_res_mat lp_map(arma::mat&) = 0;

  /* inverse of the above */
  map_res_col lp_map_inv(arma::subview_col<double> a){
   arma::vec tmp(a.colptr(0), a.n_elem, false);
     return lp_map_inv(tmp);
  }
  virtual map_res_col lp_map_inv(arma::vec&) = 0;
  virtual map_res_mat lp_map_inv(arma::mat&) = 0;
};

/* Concrete class for n-th order random walk model starting with problem data
 * for n-th order random walk */
template<class T>
class random_walk : public T {
  std::unique_ptr<arma::span> span_current_cov(){
    std::unique_ptr<arma::span> out;
    if(this->any_dynamic || this->any_fixed_in_E_step){
      out.reset(new arma::span(
          0,
          (this->space_dim - this->n_params_state_vec_fixed) / order() +
            this->n_params_state_vec_fixed - 1));
    }

    return out;
  }

  using ptr_vec = std::unique_ptr<arma::vec>;
  using ptr_mat = std::unique_ptr<arma::mat>;

public:
  using T::T;

  /* particular class function */
  unsigned int order(){
    unsigned int numerator = this->space_dim - this->n_params_state_vec_fixed;
    if(numerator == 0)
      return 1;

    return numerator / (this->covar_dim - this->n_params_state_vec_fixed);
  }

  /* derived class function to override */
  map_res_col lp_map(arma::vec &a) override {
    auto span_use = (order() == 1) ? arma::span::all : *span_current_cov();
    return map_res_col(a(span_use));
  }

  map_res_mat lp_map(arma::mat &M) override {
    auto span_use = (order() == 1) ? arma::span::all : *span_current_cov();
    return map_res_mat(M(span_use, span_use));
  }

  map_res_col lp_map_inv(arma::vec &a) override {
    if(order() == 1)
      return map_res_col(a(arma::span::all));

    ptr_vec ptr(new arma::vec(this->space_dim, arma::fill::zeros));
    arma::vec &out = *ptr.get();
    out(*span_current_cov()) = a;

    return map_res_col(out(arma::span::all), ptr);
  }

  map_res_mat lp_map_inv(arma::mat &M) override {
    if(order() == 1)
      return map_res_mat(M(arma::span::all, arma::span::all));

    ptr_mat ptr(new arma::mat(this->space_dim, this->space_dim, arma::fill::zeros));
    arma::mat &out = *ptr.get();
    out(*span_current_cov(), *span_current_cov()) = M;

    return map_res_mat(out(arma::span::all, arma::span::all), ptr);
  }
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
                const arma::mat &R,
                const arma::mat &L,
                const arma::vec &m,
                arma::mat &Q_0,
                arma::mat &Q,
                const Rcpp::List &risk_obj,
                const arma::mat &F_,
                const double eps_fixed_parems,
                const int max_it_fixed_params,
                const arma::vec &weights,
                const int n_max, const double eps,
                const bool verbose,
                const bool est_Q_0,
                const bool debug,
                Rcpp::Nullable<Rcpp::NumericVector> LR,
                const int n_threads,
                const double denom_term,
                const bool use_pinv):
    problem_data(
      n_fixed_terms_in_state_vec,
      X,
      fixed_terms,
      tstart,
      tstop, is_event_in_bin,
      a_0,
      R,
      L,
      m,
      Q_0,
      Q,
      risk_obj,
      F_,
      n_max,
      n_threads),
    weights(weights),

    event_eps(d * std::numeric_limits<double>::epsilon()),
    eps_fixed_parems(eps_fixed_parems),
    max_it_fixed_params(max_it_fixed_params),
    denom_term(denom_term),

    debug(debug),
    LR(LR.isNotNull() ? Rcpp::as< Rcpp::NumericVector >(LR)[0] : 1.0),

    use_pinv(use_pinv),
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

#endif
