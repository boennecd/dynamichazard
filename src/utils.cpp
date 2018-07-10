#include "utils.h"
#include <cmath>

/*
  We are looking for eta - exp(eta) * t = log(eps) < 0
  The answer is eta = log(eps) - W_1(-t * eps) for this function
  With eta >= -exp(eta) * t
  W_1 is defined on (-1/e, 0) yielding values in (-inf, -1) and is strictly
  decreasing
  So eta in (log(eps) + 1, inf)
  Bisection method is applied to find the solution
*/

// cache the last value
static double cache_at_risk_length =
  std::numeric_limits<double>::quiet_NaN();
static double cache_at_risk_length_result =
  std::numeric_limits<double>::quiet_NaN();

static constexpr double log_eps = trunc_eta_exponential_log_eps;

template <typename T>
inline int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

inline double func(const double eta, const double at_risk_length){
  return -(eta - exp(eta) * at_risk_length - log_eps);
}

double trunc_eta_exponential_inner_func(const double at_risk_length){
  static constexpr unsigned int n_max = 1000;
  static constexpr double delta = 10;
  static constexpr double tol = 1e-10;

  if(!std::isnan(cache_at_risk_length) &&
     cache_at_risk_length == at_risk_length)
    return cache_at_risk_length_result;

  // Find starting values for upper and lower bound for Bisection method
  double lb, ub;
  ub = log_eps + 1;
  while(func(ub, at_risk_length) < 0)
    ub += delta;

  lb = ub - delta;

  /*
    Simple Bisection method from:
      https://en.wikipedia.org/wiki/Bisection_method#Algorithm
  */
  double mb, f_mb;

  unsigned int n = 0;
  for(; n < n_max; ++n){
    mb = (lb + ub) / 2;
    f_mb = func(mb, at_risk_length);
    if(std::abs(f_mb) < tol){
      break; // solution found
    }

    if(sgn(f_mb) == sgn(func(lb, at_risk_length))){
      lb = mb;
    } else
      ub = mb;
  }

  if(n == n_max){
    std::stringstream msg;
    msg << "trunc_eta_exponential_inner_func did not converge with at_risk_length = "
        << at_risk_length;
    Rcpp::stop(msg.str());
  }

  cache_at_risk_length = at_risk_length;
  cache_at_risk_length_result = mb;

  return mb;
}

// @param x sequence of values to one of \code{boundaries} if they are almost equal.
// @param x_ord order of \code{x} in increasing order e.g. by using \code{\link{order}}. Must be a zero-based index.
// @param boundaries boundaries which \code{x} should match if it is very close to being equal.
// @details Kinda same logic as in \code{\link{all.equal.numeric}}.
// [[Rcpp::export]]
arma::vec round_if_almost_eq(
    const arma::vec &x, const arma::uvec &x_ord,
    const arma::vec &boundaries){
  const static double threshold =
    std::sqrt(std::numeric_limits<double>::epsilon());
  arma::vec x_out = x; // make copy

  double *x_out_begin = x_out.memptr();
  const double *it_bound = boundaries.begin();

  double bound_abs = std::abs(*it_bound);
  bool is_relative_test = bound_abs > threshold;
  bool is_first_it = true;
  arma::uword idx_max = x.n_elem - 1L;
  double x_old;
  for(auto i = x_ord.begin(); i != x_ord.end(); ++i){
    if(*i > idx_max)
      Rcpp::stop("`x_ord` out of bounds");

    double *this_x = x_out_begin + *i;

    /* check x_ord gives the correct order */
    if(is_first_it){
      is_first_it = false;

    } else {
      if(*this_x < x_old)
        Rcpp::stop("`x_ord` does not seem to give the correct order of `x`");

    }
    x_old = *this_x;

    /* compute test value */
    double test_val =
      is_relative_test ?
      (*this_x - *it_bound) / bound_abs :
       *this_x - *it_bound;

    if(test_val >= threshold){
      // have to proceed to next boundary and test again
      --i;
      ++it_bound;
      if(it_bound == boundaries.end())
        break;

      bound_abs = std::abs(*it_bound);
      is_relative_test = bound_abs > threshold;
      continue;
    }

    if(test_val <= -threshold) // below so no action needed
      continue;

    *this_x = *it_bound; // x is in (lb, ub) -- set x to boundary
  }

  return x_out;
}

// [[Rcpp::export]]
arma::vec rep_vec(const arma::vec &col_vals, int n_rows){
  int n_cols = col_vals.n_elem;
  arma::vec out(n_rows * n_cols);

  const double *v = col_vals.memptr();
  double *o = out.memptr();
  for(int j = 0; j < n_cols; ++j, ++v)
    for(int i = 0; i < n_rows; ++i, ++o)
      *o = *v;

  return out;
}


selection_matrix::selection_matrix(const arma::mat &A):
  n(A.n_rows), m(A.n_cols), A(A) {
  std::vector<arma::uword> idx_m_val;
  std::vector<arma::uword> idx_n_val;

  const double *a = A.memptr();
  for(arma::uword j = 0; j < m; ++j){
    bool found_one = false;

    for(arma::uword i = 0; i < n; ++i, ++a){
      if(*a < 1 - 3e-16 || 1 + 3e-16 < *a)
        continue;
      if(found_one)
        Rcpp::stop("A does not seem to be a selection matrix.");

      idx_n_val.push_back(i);
      idx_m_val.push_back(j);

      found_one = true;
    }
  }

  idx_n.reset(new arma::uvec(idx_n_val));
  idx_m.reset(new arma::uvec(idx_m_val));
}

arma::vec selection_matrix::map(const arma::vec &x) const{
  arma::vec out(n, arma::fill::zeros);

  out(*idx_n.get()) = x(*idx_m.get());

  return out;
}

arma::mat selection_matrix::map(
    const arma::mat &x, const bool is_right) const {
  if(!is_right){
    arma::mat out(n, x.n_cols, arma::fill::zeros);
    out.rows(*idx_n.get()) = x.rows(*idx_m.get());
    return out;
  }

  arma::mat out(x.n_rows, n, arma::fill::zeros);
  out.cols(*idx_n.get()) = x.cols(*idx_m.get());
  return out;
}

arma::vec selection_matrix::map_inv(const arma::vec &x) const {
  arma::vec out(m, arma::fill::zeros);
  out(*idx_m.get()) = x(*idx_n.get());

  return out;
}

arma::mat selection_matrix::map_inv(
    const arma::mat &x, const bool is_right) const {
  if(!is_right){
    arma::mat out(m, x.n_cols, arma::fill::zeros);
    out.rows(*idx_m.get()) = x.rows(*idx_n.get());
    return out;
  }

  arma::mat out(x.n_rows, m, arma::fill::zeros);
  out.cols(*idx_m.get()) = x.cols(*idx_n.get());
  return out;
}

const arma::uvec& selection_matrix::non_zero_row_idx() const {
  return *idx_n.get();
}

const arma::uvec& selection_matrix::non_zero_col_idx() const {
  return *idx_m.get();
}
