#include "problem_data.h"
#include "R_BLAS_LAPACK.h"

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
    const Rcpp::List risk_sets, const unsigned int t /* refers to bin 0, 1, 2, ... */){
  return Rcpp::as<arma::uvec>(risk_sets[t - 1]) - 1;
}

inline arma::uvec get_risk_set(
    const problem_data &data, const unsigned int t /* refers to bin 0, 1, 2, ... */){
  return get_risk_set(data.risk_sets, t);
}

struct get_bin_times_result {
  double start, stop;
};

inline get_bin_times_result get_bin_times(
    const double min_start, const std::vector<double> &I_len,
    const unsigned int t /* refers to bin 0, 1, 2, ... */){
  get_bin_times_result ans;
  double &start = ans.start;
  double &stop = ans.stop;

  start = min_start;
  unsigned int i = 0;
  for(; i + 1 < t; ++i)
    start += I_len[i];

  stop = start + I_len[i];

  return ans;
}

inline get_bin_times_result get_bin_times
  (const problem_data &data, const unsigned int t ){
  return get_bin_times(data.min_start, data.I_len, t);
}


inline double lambert_W0(const double x){
  /*
    Taylor series from:
      https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions
  */
  return x * (1 - x * (1 - 3/2 * x * (1 - 4 * x * (1 - 125 / 144 * x))));
}

struct trunc_eta_res {
  double eta_trunc;
  double exp_eta_trunc;
};

// Function to truncate the linear predictor for exponentially distributed
// outcomes
/* CHECK: The next two constants match. I have hard coded them as the exp is
   not constexpr with clang-4.0 */
static constexpr double trunc_eta_exponential_log_eps = -50;
static constexpr double trunc_eta_exponential_eps = 1.9287498479639178e-22;
double trunc_eta_exponential_inner_func(const double);

inline trunc_eta_res
  trunc_eta_exponential(
    const bool outcome, const double eta, const double exp_eta,
    const double at_risk_length)
  {
    static constexpr double log_eps = trunc_eta_exponential_log_eps;
    static constexpr double eps = trunc_eta_exponential_eps;

    trunc_eta_res ans;
    ans.exp_eta_trunc = exp_eta;

    // P(outcome) < eps or f(outcome) < eps
    bool did_truncate =
      outcome * eta - ans.exp_eta_trunc * at_risk_length < log_eps;

    if(!did_truncate){
      ans.eta_trunc = eta;
      return ans;
    }

    if(outcome){
      if(eta < -ans.exp_eta_trunc * at_risk_length){
        // answer is eta = log(eps) - W_0(-t * eps)
        ans.eta_trunc = log_eps - lambert_W0(- at_risk_length * eps);

      } else {
        // answer is eta = log(eps) - W_1(-t * eps)
        ans.eta_trunc = trunc_eta_exponential_inner_func(at_risk_length);

      }

    } else
      ans.eta_trunc = log(- log_eps / at_risk_length);

    ans.exp_eta_trunc = exp(ans.eta_trunc);
    return ans;
  }


inline arma::vec get_linear_product(
    const arma::vec &coef, const arma::mat &X, const arma::uvec &rset){
  /* Next line is slow in Armadillo release 7.960.1 (Northern Banana
     Republic Deluxe) due to memory copy in subview::extract and
     subview_elem2::extract. Same goes for similar calls */

  // arma::vec eta =  coef.t() * X.cols(rset);

  // turns out this is not much faster
  arma::vec eta(rset.n_elem);
  double *it_e = eta.begin();
  const arma::uword *it_idx = rset.begin();
  const int q = X.n_rows;
  const double *coef_ptr = coef.memptr();
  int inc = 1;
  for(unsigned int k = 0; k < rset.n_elem; ++k, ++it_e, ++it_idx)
    *it_e = R_BLAS_LAPACK::ddot(&q, coef_ptr, &inc, X.colptr(*it_idx), &inc);

  return eta;
}

#undef MIN
#undef MAX
#endif
