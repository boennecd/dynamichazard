#ifndef PF_UTILS
#define PF_UTILS

#include <tuple>
#include "particles.h"
#include "densities.h"
#include "../arma_BLAS_LAPACK.h"
#include "../utils.h"
#include "cond_approx.h"

#define MAX(a,b) (((a)>(b))?(a):(b))

struct normalize_weights_output {
  double ESS = 0.;
  double log_sum_logs;
  arma::vec weights;
};

template<class F, bool compute_ESS, bool update, typename T>
inline normalize_weights_output normalize_weights(T &container, const double max_weight){
  normalize_weights_output ans;
  double &ESS = ans.ESS;
  arma::vec &weights = ans.weights;
  weights.set_size(container.size());

  auto w = weights.begin();
  double norm_constant = 0;
  for(auto it = container.begin(); it != container.end(); ++it, ++w){
    /* back transform weights */
    *w = MAX(exp(F::get(*it) - max_weight), std::numeric_limits<double>::epsilon());

    norm_constant += *w;
  }
  ans.log_sum_logs = log(norm_constant) + max_weight;

  w = weights.begin();
  for(auto it = container.begin(); it != container.end(); ++it, ++w){
    *w /= norm_constant;

    if(compute_ESS){
      ESS += *w * *w;
    }

    if(update){
      /* Re-scale and take log */
      F::get(*it) = log(*w);
    }
  }

  if(compute_ESS){
    ESS = 1/ESS;
  }

  return ans;
}

template<class F, bool compute_ESS, bool update, typename T>
inline normalize_weights_output normalize_weights(T &container){
  double max_weight = -std::numeric_limits<double>::max();
  for(auto it = container.begin(); it != container.end(); ++it){
    max_weight = MAX(F::get(*it), max_weight);
  }

  return(normalize_weights<F, compute_ESS, update>(container, max_weight));
}

struct normalize_log_weights_F {
  template<typename T>
  static inline double& get(T &p){
    return p.log_weight;
  }
};
template
  <bool compute_ESS, bool update_particles, class TContainer>
inline normalize_weights_output normalize_log_weights(
    TContainer &container, const double max_weight){

  return
  normalize_weights
  <normalize_log_weights_F, compute_ESS, update_particles>
  (container, max_weight);
}

struct normalize_log_resampling_weight_F{
  static inline double& get(particle &p){
    return p.log_resampling_weight;
  }
};
template<bool compute_ESS, bool update_particles>
inline normalize_weights_output normalize_log_resampling_weight(
    cloud &cl, const double max_weight){
  return
    normalize_weights
    <normalize_log_resampling_weight_F, compute_ESS, update_particles>
    (cl, max_weight);
}

/* ------------------------------------------- */

struct nothing {};

/* ------------------------------------------- */

template<bool is_forward>
std::vector<std::unique_ptr<dist_comb>> get_approx_use_mean(
    std::shared_ptr<PF_cdist>, cloud&, const PF_data&, pf_dens&, arma::uword);

template<bool is_forward>
std::vector<std::unique_ptr<dist_comb>> get_approx_use_particle(
    std::shared_ptr<PF_cdist>, cloud&, const PF_data&, pf_dens&, arma::uword);

/* ------------------------------------------- */

/* Output class for smoothers */

class smoother_output {
public:
  struct pair {
    const particle *p;
    double log_weight;

    pair(const particle *p = nullptr,
         double log_weight = std::numeric_limits<double>::quiet_NaN());
  };

  struct particle_pairs {
    const particle *p;
    double log_weight;
    std::vector<pair> transition_pairs;

    particle_pairs(const particle*, const double, std::vector<pair>&&);
    particle_pairs(const particle*, const double);
    particle_pairs();
  };

  using trans_like_obj = std::vector<std::vector<particle_pairs>>;

  std::vector<cloud> forward_clouds;
  std::vector<cloud> backward_clouds;
  std::vector<cloud> smoothed_clouds;

  smoother_output();

  std::shared_ptr<trans_like_obj>
    get_transition_likelihoods(const bool do_make_if_len_0 = false) const;

private:
  std::shared_ptr<trans_like_obj> transition_likelihoods;
};

Rcpp::List get_rcpp_list_from_cloud(
    const smoother_output &sm_output, const PF_data *data = nullptr);

Rcpp::List get_rcpp_list_from_cloud(
    const std::vector<cloud> &clouds, const bool reverse,
    const unsigned int state_dim, const PF_data *data = nullptr);

smoother_output get_clouds_from_rcpp_list(const Rcpp::List &rcpp_list);

/* ------------------------------------------- */

cloud re_sample_cloud(const unsigned int, const cloud);

#undef USE_PRIOR_IN_BW_FILTER_DEFAULT
#undef MAX
#endif
