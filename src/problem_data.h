#ifndef DDHAZARD_DATA
#define DDHAZARD_DATA

#include "arma_n_rcpp.h"
#include "arma_BLAS_LAPACK.h"
#include "lin_maps.h"
#include <future>
#include <type_traits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

class problem_data {
protected:
  virtual
  std::unique_ptr<linear_mapper> set_state_trans(const arma::mat &F){
    return std::unique_ptr<dens_mapper>(new dens_mapper(F));
  }

  virtual
  std::unique_ptr<linear_mapper> set_state_trans_inv(const arma::mat &F){
    return std::unique_ptr<inv_mapper>(new inv_mapper(F));
  }

  virtual
  std::unique_ptr<linear_mapper> set_state_trans_err
  (const arma::mat &F, const arma::mat R){
    arma::mat tmp = R.t() * F;
    return std::unique_ptr<dens_mapper>(new dens_mapper(tmp));
  }

  virtual
  std::unique_ptr<linear_mapper> set_err_state(const arma::mat &R){
    return std::unique_ptr<select_mapper>(new select_mapper(R));
  }

  virtual
  std::unique_ptr<linear_mapper> set_err_state_inv(const arma::mat &R){
    /* yielded and UBSAN error when R has no columns or rows in `Rcpparmadillo`
     * version 0.8.500.0 in commit 9ae4dc1aeeba of this package
       #0 0x7f5b21117053 in arma::Mat<double>::at(unsigned int, unsigned
          int) const /home/ben/R_check_clang/lib/R/library/RcppArmadillo/include/armadillo_bits/Mat_meat.hpp:5516:3
       #1 0x7f5b21117053 in void
          arma::op_strans::apply_mat_noalias<double, arma::Mat<double>
          >(arma::Mat<double>&, arma::Mat<double> const&)
          /home/ben/R_check_clang/lib/R/library/RcppArmadillo/include/armadillo_bits/op_strans_meat.hpp:128
       #2 0x7f5b211167ac in void arma::op_strans::apply_mat<double,
          arma::Mat<double> >(arma::Mat<double>&, arma::Mat<double> const&)
          /home/ben/R_check_clang/lib/R/library/RcppArmadillo/include/armadillo_bits/op_strans_meat.hpp:212:5
       #3 0x7f5b211167ac in void
          arma::op_strans::apply_proxy<arma::Mat<double>
          >(arma::Mat<arma::Mat<double>::elem_type>&, arma::Mat<double> const&)
          /home/ben/R_check_clang/lib/R/library/RcppArmadillo/include/armadillo_bits/op_strans_meat.hpp:239
       #4 0x7f5b21199168 in void arma::op_htrans::apply<arma::Mat<double>
          >(arma::Mat<arma::Mat<double>::elem_type>&,
          arma::Op<arma::Mat<double>, arma::op_htrans> const&,
          arma::arma_not_cx<arma::Mat<double>::elem_type>::result const*)
          /home/ben/R_check_clang/lib/R/library/RcppArmadillo/include/armadillo_bits/op_htrans_meat.hpp:288:3
       #5 0x7f5b21199168 in arma::Mat<double>::Mat<arma::Mat<double>,
          arma::op_htrans>(arma::Op<arma::Mat<double>, arma::op_htrans> const&)
          /home/ben/R_check_clang/lib/R/library/RcppArmadillo/include/armadillo_bits/Mat_meat.hpp:4593
     */
    arma::mat use;
    if(R.n_cols > 0 and R.n_rows > 0)
      use = R.t();
    return std::unique_ptr<select_mapper>(new select_mapper(use));
  }

  virtual
  std::unique_ptr<linear_mapper> set_state_lp(const arma::mat &L){
    return std::unique_ptr<select_mapper>(new select_mapper(L));
  }

  virtual
  std::unique_ptr<linear_mapper> set_state_lp_inv(const arma::mat &L){
    arma::mat use;
    if(L.n_cols > 0 and L.n_rows > 0)
      use = L.t();
    return std::unique_ptr<select_mapper>(new select_mapper(arma::mat(use)));
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

  // maps
  const std::unique_ptr<linear_mapper> state_trans;
  const std::unique_ptr<linear_mapper> state_trans_inv;
  const std::unique_ptr<linear_mapper> state_trans_err;

  const std::unique_ptr<linear_mapper> err_state;
  const std::unique_ptr<linear_mapper> err_state_inv;

  /* here to allow for fixed effects in E-step */
  const std::unique_ptr<linear_mapper> state_lp;
  const std::unique_ptr<linear_mapper> state_lp_inv;

  problem_data(
    const int n_fixed_terms_in_state_vec,
    arma::mat &X,
    arma::mat &fixed_terms,
    const arma::vec &tstart,
    const arma::vec &tstop, const arma::ivec &is_event_in_bin,
    const arma::colvec &a_0,
    const arma::mat &R,
    const arma::mat &L,
    arma::mat &Q_0,
    arma::mat &Q,
    const Rcpp::List &risk_obj,
    const arma::mat &F_,
    const int n_max,
    const int n_threads,
    const arma::vec &fixed_parems) :
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
                arma::vec(X.n_cols, arma::fill::zeros)),

    state_trans(set_state_trans(F_)),
    state_trans_inv(set_state_trans_inv(F_)),
    state_trans_err(set_state_trans_err(F_, R)),

    err_state(set_err_state(R)),
    err_state_inv(set_err_state_inv(R)),

    state_lp(set_state_lp(L)),
    state_lp_inv(set_state_lp_inv(L))
  {
#ifdef _OPENMP
    omp_set_num_threads(n_threads);
    omp_set_max_active_levels(1);
#endif
  }

  problem_data& operator=(const problem_data&) = delete;
  problem_data(const problem_data&) = delete;
  problem_data() = delete;

  // create a virtual, default destructor
  virtual ~problem_data() = default;
};

/* Concrete class for n-th order random walk model starting with problem data
 * for n-th order random walk */
/* TODO: implement more overrides */

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

  std::unique_ptr<linear_mapper> set_state_trans(const arma::mat &F) override {
    if(order == 1)
      return std::unique_ptr<select_mapper>(new select_mapper(F));

    return T::set_state_trans(F);
  }

  std::unique_ptr<linear_mapper> set_state_trans_inv
  (const arma::mat &F) override {
    if(order == 1)
      return std::unique_ptr<select_mapper>(new select_mapper(F));

    return T::set_state_trans_inv(F);
  }

  std::unique_ptr<linear_mapper> set_state_trans_err
  (const arma::mat &F, const arma::mat R) override {
    if(order == 1)
      return std::unique_ptr<select_mapper>(new select_mapper(F));

    return T::set_state_trans_err(F, R);
  }

public:
  using T::T;
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
      n_fixed_terms_in_state_vec, X, fixed_terms, tstart, tstop,
      is_event_in_bin, a_0, R, L, Q_0, Q, risk_obj, F_, n_max, n_threads,
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

#endif
