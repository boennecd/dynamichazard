#include "problem_data.h"

#ifndef DDHAZARD_UTILS
#define DDHAZARD_UTILS
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

inline double get_at_risk_length(
    const double obs_stop, const double bin_stop,
    const double obs_start, const double bin_start){
  return MIN(obs_stop, bin_stop) - MAX(obs_start, bin_start);
}

inline arma::uvec get_risk_set(
    const problem_data &data, const unsigned int t /* refers to bin 0, 1, 2, ... */){
  return Rcpp::as<arma::uvec>(data.risk_sets[t - 1]) - 1;
}

struct get_bin_times_result {
  double start, stop;
};

inline get_bin_times_result get_bin_times(
    const problem_data &data, const unsigned int t /* refers to bin 0, 1, 2, ... */){
  get_bin_times_result ans;
  double &start = ans.start;
  double &stop = ans.stop;

  start = data.min_start;
  unsigned int i = 0;
  for(; i + 1 < t; ++i)
    start += data.I_len[i];

  stop = start + data.I_len[i];

  return ans;
}

inline double lambert_W0(const double x){
  /*
    Taylor series from:
      https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions
  */
  return x * (1 - x * (1 - 3/2 * x * (1 - 4 * x * (1 - 125 / 144 * x))));
}

struct trunc_lp_in_exponential_dist_result {
  double eta_trunc;
  double exp_eta_trunc;
  bool did_truncate;
};

// Function to truncate the linear predictor for exponentially distributed outcomes
static constexpr double trunc_lp_in_exponential_dist_log_eps = -150;
double trunc_lp_in_exponential_dist_inner_func(const double);

inline trunc_lp_in_exponential_dist_result
  trunc_lp_in_exponential_dist(
    const double eta, const double at_risk_length, const bool is_event)
  {
    static constexpr double log_eps = trunc_lp_in_exponential_dist_log_eps;
    static constexpr double eps = exp(log_eps);

    trunc_lp_in_exponential_dist_result ans;
    ans.exp_eta_trunc = exp(eta);

    // P(outcome) < eps or f(outcome) < eps
    ans.did_truncate = is_event * eta - ans.exp_eta_trunc * at_risk_length < log_eps;

    if(!ans.did_truncate){
      ans.eta_trunc = eta;
      return ans;
    }

    if(is_event){
      if(eta < -ans.exp_eta_trunc * at_risk_length){
        // answer is eta = log(eps) - W_0(-t * eps)
        ans.eta_trunc = log_eps - lambert_W0(- at_risk_length * eps);

      } else {
        // answer is eta = log(eps) - W_1(-t * eps)
        ans.eta_trunc = trunc_lp_in_exponential_dist_inner_func(at_risk_length);

      }

    } else
      ans.eta_trunc = log(- log_eps / at_risk_length);

    ans.exp_eta_trunc = exp(ans.eta_trunc);
    return ans;
  }

#undef MIN
#undef MAX
#endif
