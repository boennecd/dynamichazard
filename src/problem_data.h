#ifndef DDHAZARD_DATA
#define DDHAZARD_DATA

#include "arma_n_rcpp.h"
#include "arma_BLAS_LAPACK.h"
#include <future>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Object used for map function. The unique pointer need to be to a new object
 * created in the functions to make sure the orginal object is not destructed.
 */

#define PROTECTED_MAP_FUNCS(name, vec_var, map_mat)          \
virtual map_res_col name(arma::vec &a, ptr_vec &ptr) const { \
  ptr.reset(new arma::vec(map_mat * vec_var));               \
  arma::vec &out = *ptr.get();                               \
  return map_res_col(out(arma::span::all), ptr);             \
}

# define GENERIC_MAP_FUNCS_VEC(name)                        \
map_res_col name(arma::subview_col<double> a) const {       \
  ptr_vec ptr(new arma::vec(&a.at(0, 0), a.n_rows, false)); \
  return name(*ptr.get(), ptr);                             \
}                                                           \
map_res_col name(arma::vec &a) const {                      \
  ptr_vec ptr;                                              \
  return name(a, ptr);                                      \
}                                                           \
/* TODO: re-implement to avoid copy */                      \
map_res_col name(const arma::vec &a) const {                \
  ptr_vec ptr(new arma::vec(a));                            \
  return name(*ptr.get(), ptr);                             \
}

# define GENERIC_MAP_FUNCS_MAT(name, dsn_mat, transpose)            \
virtual map_res_mat name(const arma::mat &M, side s = both) const { \
  return mat_map(dsn_mat, M, s, transpose);                         \
}

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

  const LU_factorization F_LU;
  const arma::mat F_inv; // TODO: use forward and backward solves with F_LU

  using ptr_vec = std::unique_ptr<arma::vec>;
  using ptr_mat = std::unique_ptr<arma::mat>;

  /* matrix maps */
  static map_res_mat
    mat_map
    (const arma::mat &dsn_mat, const arma::mat &M, side s, bool tranpose){
      ptr_mat ptr;

      if(tranpose){
        switch(s){
        case  left:
          ptr.reset(new arma::mat(dsn_mat.t() * M));
          break;
        case both:
          ptr.reset(new arma::mat(dsn_mat.t() * M * dsn_mat));
          break;
        case right:
          ptr.reset(new arma::mat(              M * dsn_mat));
          break;
        default:
          Rcpp::stop("'Side' not implemented");
        }
      } else
        switch(s){
        case  left:
          ptr.reset(new arma::mat(dsn_mat * M));
          break;
        case both:
          ptr.reset(new arma::mat(dsn_mat * M * dsn_mat.t()));
          break;
        case right:
          ptr.reset(new arma::mat(          M * dsn_mat.t()));
          break;
        default:
          Rcpp::stop("'Side' not implemented");
        }

      arma::mat &out = *ptr.get();
      return map_res_mat(out(arma::span::all, arma::span::all), ptr);
    }

  /* vector maps */
  PROTECTED_MAP_FUNCS(lp_map, a, L)
  PROTECTED_MAP_FUNCS(lp_map_inv, a, L.t())

  PROTECTED_MAP_FUNCS(err_state_map, a, R)
  PROTECTED_MAP_FUNCS(err_state_map_inv, a, R.t())

  PROTECTED_MAP_FUNCS(state_trans_map, arma::vec(a + m), F_)
  virtual map_res_col state_trans_map_inv(arma::vec &a, ptr_vec &ptr) const {
    ptr.reset(new arma::vec(F_LU.solve(arma::vec(a - m))));
    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }

public:
  // Constants
  const bool any_dynamic;         // Any time-varying coefficients?
  const bool any_fixed_in_E_step;
  const bool any_fixed_in_M_step;

  const int d;
  const Rcpp::List risk_sets;

  const unsigned int space_dim; // dimension of state space vector
  const unsigned int covar_dim; // dimension of dynamic covariate vector
  const unsigned int err_dim;   // dimension of state space error term
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
  arma::vec fixed_parems;

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
    const int n_threads,
    const arma::vec &fixed_parems) :
    F_(F_), R(R), L(L), m(m), F_LU(F_), F_inv(arma::inv(F_)),

    any_dynamic(X.n_elem > 0),
    any_fixed_in_E_step(n_fixed_terms_in_state_vec > 0),
    any_fixed_in_M_step(fixed_terms.n_elem > 0),

    d(Rcpp::as<int>(risk_obj["d"])),
    risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"])),

    space_dim(a_0.size()),
    covar_dim(X.n_rows),
    err_dim(Q.n_cols),
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
    Q_0(Q_0),
    fixed_parems(fixed_parems)
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
  GENERIC_MAP_FUNCS_VEC(state_trans_map)
  GENERIC_MAP_FUNCS_MAT(state_trans_map, F_, false)

  /* inverse of the above */
  GENERIC_MAP_FUNCS_VEC(state_trans_map_inv)
  GENERIC_MAP_FUNCS_MAT(state_trans_map_inv, F_inv, false)

  /* maps from errors to state space dimension */
  GENERIC_MAP_FUNCS_VEC(err_state_map)
  GENERIC_MAP_FUNCS_MAT(err_state_map, R, false)

  /* inverse of the above */
  GENERIC_MAP_FUNCS_VEC(err_state_map_inv)
  GENERIC_MAP_FUNCS_MAT(err_state_map_inv, R, true)

  /* maps from state space dimension to covariate dimension */
  GENERIC_MAP_FUNCS_VEC(lp_map)
  GENERIC_MAP_FUNCS_MAT(lp_map, L, false)

  /* inverse of the above */
  GENERIC_MAP_FUNCS_VEC(lp_map_inv)
  GENERIC_MAP_FUNCS_MAT(lp_map_inv, L, true)
};

/* Concrete class for n-th order random walk model starting with problem data
 * for n-th order random walk */
/* TODO: implement more overrides */

# define RANDOM_WALK_VEC_MAP_TO_STATE(name, span_func)                       \
map_res_col name(arma::vec &a, ptr_vec &ptr) const override {                \
  if(order() == 1 && !this->any_fixed_in_E_step)                             \
    return map_res_col(a(arma::span::all), ptr);                             \
                                                                             \
  ptr_vec ptr_new(new arma::vec(this->space_dim, arma::fill::zeros));        \
  arma::vec &out = *ptr_new.get();                                           \
  auto span_use_ptr = span_func();                                           \
  if(span_use_ptr)                                                           \
    out(*span_use_ptr) = a;                                                  \
                                                                             \
  return map_res_col(out(arma::span::all), ptr_new);                         \
}

# define RANDOM_WALK_VEC_MAP_FROM_STATE(name, span_func)                     \
map_res_col name(arma::vec &a, ptr_vec &ptr) const override {                \
  if(order() == 1 && !this->any_fixed_in_E_step)                             \
    return map_res_col(a(arma::span::all), ptr);                             \
                                                                             \
  auto span_use_ptr = span_func();                                           \
  if(!span_use_ptr){                                                         \
    ptr_vec ptr_new(new arma::vec());                                        \
    return map_res_col((*ptr_new)(arma::span::all), ptr_new);                \
  }                                                                          \
                                                                             \
  return map_res_col(a(*span_use_ptr), ptr);                                 \
}

# define RANDOM_WALK_MAT_MAP_FROM_STATE(name, span_func)                   \
map_res_mat name(const arma::mat &M, side s = both) const override {       \
  if(s != both)                                                            \
    Rcpp::stop("'Side' not implemented");                                  \
                                                                           \
  if(order() == 1 && !this->any_fixed_in_E_step)                           \
    return map_res_mat(M(arma::span::all, arma::span::all));               \
                                                                           \
  auto span_use_ptr = span_func();                                         \
  if(!span_use_ptr){                                                       \
    ptr_mat ptr_new(new arma::mat());                                      \
    return map_res_mat(                                                    \
        (*ptr_new)(arma::span::all, arma::span::all), ptr_new);            \
  }                                                                        \
                                                                           \
  return map_res_mat(M(*span_use_ptr, *span_use_ptr));                     \
}

# define RANDOM_WALK_MAT_MAP_TO_STATE(name, span_func)                     \
map_res_mat name(const arma::mat &M, side s = both) const override {       \
  if(s != both)                                                            \
    Rcpp::stop("'Side' not implemented");                                  \
  if(order() == 1 && !this->any_fixed_in_E_step)                           \
    return map_res_mat(M(arma::span::all, arma::span::all));               \
                                                                           \
  ptr_mat ptr(                                                             \
      new arma::mat(this->space_dim, this->space_dim, arma::fill::zeros)); \
  arma::mat &out = *ptr.get();                                             \
  auto span_use_ptr = span_func();                                         \
  if(span_use_ptr)                                                         \
    out(*span_use_ptr, *span_use_ptr) = M;                                 \
                                                                           \
  return map_res_mat(out(arma::span::all, arma::span::all), ptr);          \
}

template<class T>
class random_walk : public T {
  /* Use std::unique_ptr as Armadillo does not have an "empty" span. A null
   * pointer means empty.
   * TODO: make private object once and avoid function calls                */
  std::unique_ptr<arma::span> span_current_cov() const {
    std::unique_ptr<arma::span> out;
    int b = (this->space_dim - this->n_params_state_vec_fixed) / order() +
      this->n_params_state_vec_fixed - 1;
    if(b >= 0 && (this->any_dynamic || this->any_fixed_in_E_step))
      out.reset(new arma::span(0, b));

    return out;
  }

  std::unique_ptr<arma::span> span_current_cov_varying_only() const {
    std::unique_ptr<arma::span> out;
    int b = (this->space_dim - this->n_params_state_vec_fixed) / order() - 1;
    if(b >= 0 && this->any_dynamic)
      out.reset(new arma::span(0, b));

    return out;
  }

  using ptr_vec = std::unique_ptr<arma::vec>;
  using ptr_mat = std::unique_ptr<arma::mat>;

  /* derived class functions to override */
  RANDOM_WALK_VEC_MAP_TO_STATE(err_state_map, span_current_cov_varying_only)
  RANDOM_WALK_VEC_MAP_FROM_STATE(err_state_map_inv, span_current_cov_varying_only)
  RANDOM_WALK_VEC_MAP_FROM_STATE(lp_map, span_current_cov)
  RANDOM_WALK_VEC_MAP_TO_STATE(lp_map_inv, span_current_cov)

public:
  using T::T;

  /* particular class function */
  unsigned int order() const {
    unsigned int numerator = this->space_dim - this->n_params_state_vec_fixed;
    if(numerator == 0)
      return 1;

    return numerator / (this->covar_dim - this->n_params_state_vec_fixed);
  }

  /* derived class member functions to override */
  RANDOM_WALK_MAT_MAP_TO_STATE(err_state_map, span_current_cov_varying_only)
  RANDOM_WALK_MAT_MAP_FROM_STATE(err_state_map_inv, span_current_cov_varying_only)
  RANDOM_WALK_MAT_MAP_FROM_STATE(lp_map, span_current_cov)
  RANDOM_WALK_MAT_MAP_TO_STATE(lp_map_inv, span_current_cov)
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
      n_threads,
      fixed_parems_start),
    weights(weights),

    event_eps(d * std::numeric_limits<double>::epsilon()),
    eps_fixed_parems(eps_fixed_parems),
    max_it_fixed_params(max_it_fixed_params),
    denom_term(denom_term),

    debug(debug),
    LR(LR.isNotNull() ? Rcpp::as< Rcpp::NumericVector >(LR)[0] : 1.0),

    use_pinv(use_pinv)
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

#undef PROTECTED_MAP_FUNCS
#undef GENERIC_MAP_FUNCS_VEC
#undef GENERIC_MAP_FUNCS_MAT
#undef RANDOM_WALK_VEC_MAP_TO_STATE
#undef RANDOM_WALK_VEC_MAP_FROM_STATE
#undef RANDOM_WALK_MAT_MAP_TO_STATE
#undef RANDOM_WALK_MAT_MAP_FROM_STATE
#endif
