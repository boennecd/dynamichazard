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

static constexpr double log_eps = trunc_lp_in_exponential_dist_log_eps;

template <typename T>
inline int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

inline double func(const double eta, const double at_risk_length){
  return -(eta - exp(eta) * at_risk_length - log_eps);
}

double trunc_lp_in_exponential_dist_inner_func(const double at_risk_length){
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
    msg << "trunc_lp_in_exponential_dist_inner_func did not converge with at_risk_length = "
        << at_risk_length;
    Rcpp::stop(msg.str());
  }

  cache_at_risk_length = at_risk_length;
  cache_at_risk_length_result = mb;

  return mb;
}
