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

#undef MIN
#undef MAX
#endif
