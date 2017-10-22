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
  for(auto i = x_ord.begin(); i != x_ord.end(); ++i){
    double *this_x = x_out_begin + *i;
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
