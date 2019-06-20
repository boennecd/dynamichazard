#ifndef PF_UTILS
#define PF_UTILS

#include "particles.h"
#include "densities.h"
#include "../utils.h"
#include "cond_approx.h"
#include <set>

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
    *w = exp(F::get(*it) - max_weight);

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

struct normalize_log_resampling_weight_F {
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

struct normalize_log_std_vec_double {
  static inline double& get(double &d){
    return d;
  }
};
template<bool compute_ESS, bool update>
inline normalize_weights_output normalize_log_weights(
    std::vector<double> &ws, const double max_weight){
  return
  normalize_weights
  <normalize_log_std_vec_double, compute_ESS, update>
  (ws, max_weight);
}

/* ------------------------------------------- */

struct nothing {};

/* ------------------------------------------- */

struct get_approx_use_mean_output {
  std::vector<std::unique_ptr<dist_comb>> dists;
};

template<bool is_forward>
get_approx_use_mean_output get_approx_use_mean(
    std::shared_ptr<PF_cdist>, cloud&, const PF_data&, pf_dens&, arma::uword);

struct get_approx_use_particle_output {
  std::vector<std::unique_ptr<dist_comb>> dists;
};
template<bool is_forward>
get_approx_use_particle_output get_approx_use_particle(
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

template<bool is_smooth, const bool reverse>
std::vector<cloud> get_cloud_from_rcpp_list
  (const Rcpp::List &, const std::vector<cloud> *fw = nullptr,
   const std::vector<cloud> *bw = nullptr);

/* ------------------------------------------- */

cloud re_sample_cloud(const unsigned int, const cloud);

/* ------------------------------------------- */

/* function used when sub-sampling has been performed. Returns a map with
 * indicies of those that are still included and their adjusted weight */
template
  <class TContainer,
   double (*FWeight)(const typename TContainer::value_type&),
   double (*FReWeight)(const typename TContainer::value_type&)>
  std::map<arma::uword, double>
  get_resample_idx_n_log_weight
    (const TContainer &container,
     const arma::uvec &resample_idx)
  {
    std::map<arma::uword, double> out;
    std::map<arma::uword, unsigned int> count;
    for(auto x : resample_idx){
      std::map<arma::uword, double>::iterator it = out.find(x);
      if(it == out.end()){
        out[x] = FWeight(container[x]) - FReWeight(container[x]);
        count[x] = 1L;
        continue;

      }

      count[x] += 1L;

    }

    double max_weight = -std::numeric_limits<double>::infinity();
    for(auto &x : count){
      out[x.first] += std::log(x.second);
      max_weight = std::max(max_weight, out[x.first]);
    }

    /* re-normalize */
    double norm_constant = 0;
    for(auto &x : out){
      x.second = exp(x.second - max_weight);

      norm_constant += x.second;
    }

    const double log_norm_constant = log(norm_constant);
    for(auto &x : out)
      x.second = log(x.second) - log_norm_constant;

    return out;
  }

double get_weight_from_particle(const particle&);
double get_resample_weight_from_particle(const particle&);

/* ------------------------------------------- */

std::vector<std::set<arma::uword> > get_ancestors
  (const std::vector<cloud>&);

class score_n_hess_base {
public:
  virtual const arma::vec &get_score() const = 0;
  /* the betas in https://doi.org/10.1093/biomet/asq062 */
  virtual const arma::mat &get_hess_terms() const = 0;
  virtual double get_weight() const = 0;

  virtual ~score_n_hess_base() = default;
};

std::vector<std::unique_ptr<score_n_hess_base> > PF_get_score_n_hess
  (const std::vector<cloud>&, const arma::mat&,const arma::mat&,
   const std::vector<arma::uvec>&, const arma::ivec&, const arma::vec&,
   const arma::mat&, const arma::mat&, const arma::vec&, const arma::vec&,
   const arma::vec&, const std::string, const int, const bool, const bool);

std::vector<std::unique_ptr<score_n_hess_base> > PF_get_score_n_hess_O_N_sq
  (arma::mat&, const arma::mat&, const std::vector<arma::uvec>&,
   const Rcpp::List&, const arma::ivec&, const arma::vec&,
   arma::mat&, arma::mat&, const arma::vec&,
   const arma::vec&, const arma::vec&, const std::string, const int,
   const bool, const bool, const arma::vec&, const arma::mat&, arma::mat&,
   const arma::mat&, const arma::uword, const arma::uword, const double,
   const double, const double, Rcpp::Nullable<Rcpp::NumericVector>,
   const std::string);

#undef USE_PRIOR_IN_BW_FILTER_DEFAULT
#undef MAX
#endif
