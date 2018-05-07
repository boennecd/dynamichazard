#ifndef DDHAZARD_DATA
#define DDHAZARD_DATA

#include "arma_n_rcpp.h"
#include "arma_BLAS_LAPACK.h"
#include "lin_maps.h"
#include <future>
#include <type_traits>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Object used for map function. The unique pointer need to be to a new object
 * created in the functions to make sure the orginal object is not destructed.
 */

#define PROTECTED_MAP_FUNCS(name, vec_var, Smat)                              \
virtual map_res_col name(arma::vec &a, ptr_vec &ptr) const {                  \
  ptr.reset(new arma::vec(Smat.map(vec_var)));                                \
  arma::vec &out = *ptr.get();                                                \
  return map_res_col(out(arma::span::all), ptr);                              \
}                                                                             \
virtual map_res_col name##_inv(arma::vec &a, ptr_vec &ptr) const {            \
  ptr.reset(new arma::vec(Smat.map_inv(vec_var)));                            \
  arma::vec &out = *ptr.get();                                                \
  return map_res_col(out(arma::span::all), ptr);                              \
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

/* problem_data is data holder class. It has some virtual map functions which
 * can be overwritten for a particular class of problems. The methods are
 * lp_map and lp_map_inv
 *   returns e.g., (L alpha, L A L^\top) and (L^\top x, L^\top X L).
 * err_state_map and err_state_map
 *   returns e.g., (R epsilon, R Q R^\top) and (R^-1 alpha, R^\top A R)
 * state_trans_map and state_trans_map_inv
 *   returns e.g., (F alpha + m, F A F^\top) and (F^-1(alpha - m),
 *   F^-1 A (F^top)^-1)
 *
 * It is assumed that L and R are potentially non-square matrices consisting of
 * unit vectors. Further, F is assumed to be a non-singular square matrix.
 *
 * TODO: deal with more general L and R matrices or exploit the structure to
 *       to reduce the computation time                                      */

class problem_data {
protected:
  const arma::mat &F_;
  const arma::mat &R;
  const arma::mat &L;
  const arma::vec &m; // TODO: remove

  const selection_matrix R_Smat;
  const selection_matrix L_Smat;
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
  static map_res_mat
    mat_map
    (const selection_matrix &Smat, const arma::mat &M, side s, bool is_inv){
      ptr_mat ptr;

      if(is_inv){
        switch(s){
        case  left:
          ptr.reset(new arma::mat(             Smat.map_inv(M)       ));
          break;
        case both:
          ptr.reset(new arma::mat(Smat.map_inv(Smat.map_inv(M), true)));
          break;
        case right:
          ptr.reset(new arma::mat(Smat.map_inv(             M , true)));
          break;
        default:
          Rcpp::stop("'Side' not implemented");
        }
      } else
        switch(s){
        case  left:
          ptr.reset(new arma::mat(         Smat.map(M)       ));
          break;
        case both:
          ptr.reset(new arma::mat(Smat.map(Smat.map(M), true)));
          break;
        case right:
          ptr.reset(new arma::mat(Smat.map(         M , true)));
          break;
        default:
          Rcpp::stop("'Side' not implemented");
        }

      arma::mat &out = *ptr.get();
      return map_res_mat(out(arma::span::all, arma::span::all), ptr);
    }

  /* vector maps */
  PROTECTED_MAP_FUNCS(lp_map, a, L_Smat)
  PROTECTED_MAP_FUNCS(err_state_map, a, R_Smat)

  virtual map_res_col state_trans_map(arma::vec &a, ptr_vec &ptr) const {
    ptr.reset(new arma::vec(F_ * a + m));
    arma::vec &out = *ptr.get();
    return map_res_col(out(arma::span::all), ptr);
  }
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

  const unsigned int state_dim; // dimension of state space vector
  const unsigned int covar_dim; // dimension of dynamic covariate vector
  const unsigned int err_dim;   // dimension of state space error term
  const int n_params_state_vec_fixed; // TODO: maybe move to derived class

  /* these are not const due the arma::mat constructor which does not allow
   * copy_aux_mem = false with const pointer. They should not be changed
   * though...                                                              */
  arma::mat X;
  arma::mat fixed_terms;

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
  arma::vec fixed_effects;

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
    F_(F_), R(R), L(L), m(m), R_Smat(R), L_Smat(L), F_LU(F_),
    F_inv(arma::inv(F_)),

    any_dynamic(X.n_elem > 0),
    any_fixed_in_E_step(n_fixed_terms_in_state_vec > 0),
    any_fixed_in_M_step(fixed_terms.n_elem > 0),

    d(Rcpp::as<int>(risk_obj["d"])),
    risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"])),

    state_dim(a_0.size()),
    covar_dim(X.n_rows),
    err_dim(Q.n_cols),
    n_params_state_vec_fixed(n_fixed_terms_in_state_vec),

    X(X.begin(), X.n_rows, X.n_cols, false),
    fixed_terms(fixed_terms.begin(), fixed_terms.n_rows,
                fixed_terms.n_cols, false),
    I_len(Rcpp::as<std::vector<double> >(risk_obj["I_len"])),

    n_threads((n_threads > 0) ? n_threads : std::thread::hardware_concurrency()),

    tstart(tstart),
    tstop(tstop),
    is_event_in_bin(is_event_in_bin),
    min_start(Rcpp::as<double>(risk_obj["min_start"])),

    Q(Q),
    Q_0(Q_0),
    fixed_parems(fixed_parems),
    fixed_effects(
      (any_fixed_in_M_step) ?
                fixed_terms.t() * fixed_parems :
                arma::vec(X.n_cols, arma::fill::zeros))
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
  GENERIC_MAP_FUNCS_MAT(err_state_map, R_Smat, false)

  /* inverse of the above */
  GENERIC_MAP_FUNCS_VEC(err_state_map_inv)
  GENERIC_MAP_FUNCS_MAT(err_state_map_inv, R_Smat, true)

  /* maps from state space dimension to covariate dimension */
  GENERIC_MAP_FUNCS_VEC(lp_map)
  GENERIC_MAP_FUNCS_MAT(lp_map, L_Smat, false)

  /* inverse of the above */
  GENERIC_MAP_FUNCS_VEC(lp_map_inv)
  GENERIC_MAP_FUNCS_MAT(lp_map_inv, L_Smat, true)
};

/* Concrete class for n-th order random walk model starting with problem data
 * for n-th order random walk */
/* TODO: implement more overrides */

# define RANDOM_WALK_VEC_MAP_TO_STATE(name, span_use_ptr)                    \
map_res_col name(arma::vec &a, ptr_vec &ptr) const override {                \
  if(order == 1 && !this->any_fixed_in_E_step)                               \
    return map_res_col(a(arma::span::all), ptr);                             \
                                                                             \
  ptr_vec ptr_new(new arma::vec(this->state_dim, arma::fill::zeros));        \
  arma::vec &out = *ptr_new.get();                                           \
  if(span_use_ptr)                                                           \
    out(*span_use_ptr) = a;                                                  \
                                                                             \
  return map_res_col(out(arma::span::all), ptr_new);                         \
}

# define RANDOM_WALK_VEC_MAP_FROM_STATE(name, span_use_ptr)                  \
map_res_col name(arma::vec &a, ptr_vec &ptr) const override {                \
  if(order == 1 && !this->any_fixed_in_E_step)                               \
    return map_res_col(a(arma::span::all), ptr);                             \
                                                                             \
  if(!span_use_ptr){                                                         \
    ptr_vec ptr_new(new arma::vec());                                        \
    return map_res_col((*ptr_new)(arma::span::all), ptr_new);                \
  }                                                                          \
                                                                             \
  return map_res_col(a(*span_use_ptr), ptr);                                 \
}

# define RANDOM_WALK_MAT_MAP_FROM_STATE(name, span_use_ptr)                \
map_res_mat name(const arma::mat &M, side s = both) const override {       \
  if(s != both)                                                            \
    Rcpp::stop("'Side' not implemented");                                  \
                                                                           \
  if(order == 1 && !this->any_fixed_in_E_step)                             \
    return map_res_mat(M(arma::span::all, arma::span::all));               \
                                                                           \
  if(!span_use_ptr){                                                       \
    ptr_mat ptr_new(new arma::mat());                                      \
    return map_res_mat(                                                    \
      (*ptr_new)(arma::span::all, arma::span::all), ptr_new);              \
  }                                                                        \
                                                                           \
  return map_res_mat(M(*span_use_ptr, *span_use_ptr));                     \
}

# define RANDOM_WALK_MAT_MAP_TO_STATE(name, span_use_ptr)                  \
map_res_mat name(const arma::mat &M, side s = both) const override {       \
  if(s != both)                                                            \
    Rcpp::stop("'Side' not implemented");                                  \
  if(order == 1 && !this->any_fixed_in_E_step)                             \
    return map_res_mat(M(arma::span::all, arma::span::all));               \
                                                                           \
  ptr_mat ptr(                                                             \
      new arma::mat(this->state_dim, this->state_dim, arma::fill::zeros)); \
  arma::mat &out = *ptr.get();                                             \
  if(span_use_ptr)                                                         \
    out(*span_use_ptr, *span_use_ptr) = M;                                 \
                                                                           \
  return map_res_mat(out(arma::span::all, arma::span::all), ptr);          \
}

template<class T>
class random_walk : public T {
public:
  const unsigned int order {
    (this->state_dim - this->n_params_state_vec_fixed == 0) ?
    1 :
    (this->state_dim - this->n_params_state_vec_fixed) /
      (this->covar_dim - this->n_params_state_vec_fixed)
  };

private:
  /* Use std::unique_ptr as Armadillo does not have an "empty" span. A null
   * pointer means empty.                                                   */
  const std::unique_ptr<arma::span> span_current_cov {
    get_span_current_cov() };
  std::unique_ptr<arma::span> get_span_current_cov() const {
    std::unique_ptr<arma::span> out;
    int b = (this->state_dim - this->n_params_state_vec_fixed) / order +
      this->n_params_state_vec_fixed - 1;
    if(b >= 0 && (this->any_dynamic || this->any_fixed_in_E_step))
      out.reset(new arma::span(0, b));

    return out;
  }

  const std::unique_ptr<arma::span> span_current_cov_varying_only {
    get_span_current_cov_varying_only() };
  std::unique_ptr<arma::span> get_span_current_cov_varying_only() const {
    std::unique_ptr<arma::span> out;
    int b = (this->state_dim - this->n_params_state_vec_fixed) / order - 1;
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

  map_res_col state_trans_map(arma::vec &a, ptr_vec &ptr) const override {
    if(order == 1)
      return map_res_col(a(arma::span::all), ptr);

    return T::state_trans_map(a, ptr);
  }
  map_res_col state_trans_map_inv(arma::vec &a, ptr_vec &ptr) const override {
    if(order == 1)
      return map_res_col(a(arma::span::all), ptr);

    return T::state_trans_map_inv(a, ptr);
  }

public:
  using T::T;

  /* derived class member functions to override */
  RANDOM_WALK_MAT_MAP_TO_STATE(err_state_map, span_current_cov_varying_only)
  RANDOM_WALK_MAT_MAP_FROM_STATE(err_state_map_inv, span_current_cov_varying_only)
  RANDOM_WALK_MAT_MAP_FROM_STATE(lp_map, span_current_cov)
  RANDOM_WALK_MAT_MAP_TO_STATE(lp_map_inv, span_current_cov)

  map_res_mat state_trans_map(
    const arma::mat &M, side s = both) const override {
    if(order == 1)
      return map_res_mat(M(arma::span::all, arma::span::all));

    return T::state_trans_map(M, s);
  }
  map_res_mat state_trans_map_inv(
    const arma::mat &M, side s = both) const override {
    if(order == 1)
      return map_res_mat(M(arma::span::all, arma::span::all));

    return T::state_trans_map_inv(M, s);
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
      a_t_t_s = arma::mat(state_dim, d + 1);
      a_t_less_s = arma::mat(state_dim, d);
      V_t_t_s = arma::cube(state_dim, state_dim, d + 1);
      V_t_less_s = arma::cube(state_dim, state_dim, d);
      B_s = arma::cube(state_dim, state_dim, d);

      a_t_t_s.col(0) = a_0;

      lag_one_cov = arma::cube(state_dim, state_dim, d);
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
