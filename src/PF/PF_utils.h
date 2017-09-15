#ifndef PF_UTILS
#define PF_UTILS

#include <tuple>
#include "particles.h"
#include "../arma_BLAS_LAPACK.h"
#include "../utils.h"

#define MAX(a,b) (((a)>(b))?(a):(b))
#define USE_PRIOR_IN_BW_FILTER_DEFAULT true

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
  return normalize_weights<normalize_log_resampling_weight_F, compute_ESS, update_particles>(cl, max_weight);
}

/* ------------------------------------------- */

struct nothing {};

/* ------------------------------------------- */

template<typename iterator>
class work_block {
public:
  iterator start;
  iterator end;
  const unsigned long block_size;

  work_block(iterator start, iterator end, const unsigned long block_size):
    start(start), end(end), block_size(block_size) {}
};

template<typename iterator>
std::vector<work_block<iterator>> get_work_blocks(
    iterator begin, iterator end, const unsigned long block_size){
  unsigned long const length = std::distance(begin, end);
  unsigned long const num_blocks= (length + block_size - 1) / block_size;

  std::vector<work_block<iterator>> ans;
  ans.reserve(num_blocks);

  iterator block_start = begin;
  for(unsigned long i = 0; i < num_blocks - 1; ++i){
    iterator block_end = block_start;
    std::advance(block_end, block_size);

    ans.emplace_back(block_start, block_end, block_size);
    block_start = block_end;
  }

  ans.emplace_back(block_start, end, std::distance(block_start, end));

  return ans;
}

/* ------------------------------------------- */

struct input_for_normal_apprx {
  arma::vec mu;
  arma::mat Sigma_inv_chol;
  arma::mat Sigma_chol;
  arma::mat sigma_chol_inv;
};

template<typename densities, unsigned int debug_lvl, bool multithread>
static input_for_normal_apprx compute_mu_n_Sigma_from_normal_apprx(
    const PF_data &data,
    const unsigned int t,
    const covarmat &Q,
    const arma::vec &alpha_bar){
  /*
    Had similar issues as posted here:
      http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-June/005968.html

    Thus, I made this overload
  */

  arma::uvec r_set = get_risk_set(data, t);

  return(
    compute_mu_n_Sigma_from_normal_apprx
    <densities, debug_lvl, multithread>
    (data, t, Q, alpha_bar, r_set));
}

template<typename densities, unsigned int debug_lvl, bool multithread>
static input_for_normal_apprx compute_mu_n_Sigma_from_normal_apprx(
    const PF_data &data,
    const unsigned int t,
    const covarmat &Q,
    const arma::vec &alpha_bar,
    arma::uvec &r_set){
  if(data.debug > debug_lvl){
    data.log(debug_lvl) << "Computing normal approximation with mean vector:" << std::endl
                        << alpha_bar.t()
                        << "and chol(covariance):"  << std::endl
                        << Q.chol;
  }

  input_for_normal_apprx ans;

  double bin_start, bin_stop;
  if(densities::uses_at_risk_length){
    auto tmp = get_bin_times(data, t);
    bin_start = tmp.start;
    bin_stop = tmp.stop;

  }

  /* Compute the terms that does not depend on the outcome */
  /* Sigma^-1 = (Q + \tilde{Q})^{-1} */
  arma::uword p = alpha_bar.n_elem;
  arma::mat Sigma_inv = Q.inv;

  ans.mu = arma::vec(p, arma::fill::zeros);
  arma::vec &mu = ans.mu;

  /* Add the terms that does depend on the outcome */
  auto jobs = get_work_blocks(r_set.begin(), r_set.end(), data.work_block_size);
  unsigned int n_jobs = jobs.size();

#ifdef _OPENMP
  /*
    Use lock as critical section will not do if this function is called in
    nested parallel setup. See https://stackoverflow.com/a/20447843
  */
  omp_lock_t *lock;
  if(multithread){
    lock = new omp_lock_t;
    omp_init_lock(lock);
  }
#pragma omp parallel if(multithread)
{
#endif

  arma::vec starts;
  arma::vec stops;
  arma::mat my_Sigma_inv(p, p, arma::fill::zeros);
  arma::vec my_mu(p, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for(unsigned int i = 0; i < n_jobs; ++i){
    auto &job = jobs[i];
    arma::uvec my_r_set(job.start, job.block_size, false /* don't copy */);

    arma::vec eta =  get_linear_product(alpha_bar, data.X, my_r_set);
    const arma::uvec is_event = data.is_event_in_bin(my_r_set) == t - 1; /* zero indexed while t is not */

    if(densities::uses_at_risk_length){
      starts = data.tstart(my_r_set);
      stops = data.tstop(my_r_set);
    }

    auto it_eta = eta.begin();
    auto it_is_event = is_event.begin();
    auto it_r = my_r_set.begin();
    auto it_start = starts.begin();
    auto it_stops = stops.begin();
    arma::uword n_elem = eta.n_elem;
    /*
      Update with:
        Signa = ... + X^T (-G) X
        mu = X^T (-G) X \bar{alpha} + X^T (-g)
    */
    for(arma::uword i = 0; i < n_elem; ++i, ++it_eta, ++it_is_event, ++it_r){
      double at_risk_length = 0;
      if(densities::uses_at_risk_length){
        at_risk_length = get_at_risk_length(
          *(it_stops++) /* increament here */, bin_stop,
          *(it_start++) /* increament here */, bin_start);

      }

      double g = densities::log_p_prime(*it_eta, *it_is_event, at_risk_length);
      double neg_G = - densities::log_p_2prime(*it_eta, *it_is_event, at_risk_length);

      sym_mat_rank_one_update(neg_G, data.X.col(*it_r), my_Sigma_inv);

      my_mu += data.X.col(*it_r) * ((*it_eta * neg_G) + g);
    }
  }

#ifdef _OPENMP
  if(multithread)
    omp_set_lock(lock);
#endif

  Sigma_inv += my_Sigma_inv;
  mu += my_mu;

#ifdef _OPENMP
  if(multithread)
    omp_unset_lock(lock);
#endif

#ifdef _OPENMP
}
  if(multithread){
    omp_destroy_lock(lock);
    delete lock;
  }
#endif

  /* copy to lower */
  Sigma_inv = arma::symmatu(Sigma_inv);

  /* Compute needed factorizations */
  ans.Sigma_inv_chol = arma::chol(Sigma_inv);
  ans.Sigma_chol = arma::chol(arma::inv(Sigma_inv)); // TODO: do something smarter
  ans.sigma_chol_inv = arma::inv(arma::trimatu(ans.Sigma_chol));

  return ans;
}

/* ------------------------------------------- */

struct input_for_normal_apprx_w_cloud_mean : public input_for_normal_apprx {
  /* the conditional means for each of the particles given their parents */
  std::vector<arma::vec> mu_js;

  input_for_normal_apprx_w_cloud_mean(input_for_normal_apprx &&other):
    input_for_normal_apprx(other) {}
};

template<typename densities, bool is_forward, bool use_prior_in_bw_filter = USE_PRIOR_IN_BW_FILTER_DEFAULT>
static input_for_normal_apprx_w_cloud_mean
  compute_mu_n_Sigma_from_normal_apprx_w_cloud_mean(
    const PF_data &data, const unsigned int t, const covarmat &Q, const arma::vec &alpha_bar,
    cloud &cl /* set mu_js when cloud is passed to */){
    constexpr bool is_bw_w_use_prior = (!is_forward) && use_prior_in_bw_filter;

    const covarmat *Q_use;
    const arma::mat tmp;
    arma::vec *mu_term;
    if(is_bw_w_use_prior){
      // Add the covariance matrix of the artificial prior
      arma::mat art = densities::get_artificial_prior_covar(data, t);
      Q_use = new covarmat(arma::inv(Q.inv + arma::inv(art)));
      mu_term = new arma::vec(arma::solve(art, data.a_0));

    } else {
      Q_use = &Q;

    }

    input_for_normal_apprx_w_cloud_mean ans =
      compute_mu_n_Sigma_from_normal_apprx
      <densities, 2, true>
      (data, t, *Q_use, alpha_bar);

    auto n_elem = cl.size();
    ans.mu_js = std::vector<arma::vec>(n_elem);
#ifdef _OPENMP
#pragma omp  parallel for schedule(static)
#endif
    for(unsigned int i = 0; i < n_elem; ++i){
      arma::vec mu_j = solve_w_precomputed_chol(Q.chol, cl[i].get_state()) + ans.mu;
      if(is_bw_w_use_prior)
        mu_j += *mu_term;
      mu_j = solve_w_precomputed_chol(ans.Sigma_inv_chol, mu_j);

      ans.mu_js[i] = std::move(mu_j);
    }

    if(is_bw_w_use_prior){
      delete Q_use;
      delete mu_term;

    }

    return ans;
}

/* ------------------------------------------- */

struct input_for_normal_apprx_w_particle_mean_element {
  arma::vec mu;
  arma::mat sigma_chol_inv;
  arma::mat Sigma_chol;
};

using input_for_normal_apprx_w_particle_mean =
  std::vector<input_for_normal_apprx_w_particle_mean_element>;

template<typename densities, typename mu_iterator,
         typename Func, bool is_forward, bool use_prior_in_bw_filter = USE_PRIOR_IN_BW_FILTER_DEFAULT>
static input_for_normal_apprx_w_particle_mean
compute_mu_n_Sigma_from_normal_apprx_w_particles(
  const PF_data &data, const unsigned int t, const covarmat &Q,
  mu_iterator begin, const unsigned int size){
  constexpr bool is_bw_w_use_prior = (!is_forward) && use_prior_in_bw_filter;

  const covarmat *Q_use;
  arma::mat Q_art_chol;
  arma::vec *mu_term;
  if(is_bw_w_use_prior){
    // Add the covariance matrix of the artificial prior
    arma::mat art = densities::get_artificial_prior_covar(data, t);
    Q_use = new covarmat(arma::inv(Q.inv + arma::inv(art)));
    mu_term = new arma::vec(arma::solve(art, data.a_0));

  } else {
    Q_use = &Q;

  }

  input_for_normal_apprx_w_particle_mean ans(size);

  mu_iterator b = begin;
  arma::uvec r_set = get_risk_set(data, t);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned int i = 0; i < size; ++i){
    mu_iterator iter = b + i;
    arma::vec mu = Func::get_elem(iter);
    auto inter = compute_mu_n_Sigma_from_normal_apprx
      <densities, 5, false>(data, t, *Q_use, mu, r_set);

    mu = solve_w_precomputed_chol(Q.chol, mu) + inter.mu;
    if(is_bw_w_use_prior)
      mu += *mu_term;

    mu = solve_w_precomputed_chol(inter.Sigma_inv_chol, mu);

    ans[i].mu = std::move(mu);
    ans[i].sigma_chol_inv = std::move(inter.sigma_chol_inv);
    ans[i].Sigma_chol = std::move(inter.Sigma_chol);
  }

  if(is_bw_w_use_prior){
    delete Q_use;
    delete mu_term;

  }

  return ans;
}

template<typename densities, bool is_forward, bool use_prior_in_bw_filter = USE_PRIOR_IN_BW_FILTER_DEFAULT>
static input_for_normal_apprx_w_particle_mean
compute_mu_n_Sigma_from_normal_apprx_w_particles(
  const PF_data &data, const unsigned int t, const covarmat &Q,
  cloud &cl){
  struct Func{
    static inline const arma::vec get_elem(cloud::iterator &it){
      return it->get_state();
    }
  };

  return(
    compute_mu_n_Sigma_from_normal_apprx_w_particles
    <densities, cloud::iterator, Func, is_forward, use_prior_in_bw_filter>
    (data, t, Q, cl.begin(), cl.size()));
}

template<typename densities, bool is_forward, bool use_prior_in_bw_filter = USE_PRIOR_IN_BW_FILTER_DEFAULT>
static input_for_normal_apprx_w_particle_mean
compute_mu_n_Sigma_from_normal_apprx_w_particles(
  const PF_data &data, const unsigned int t, const covarmat &Q,
  std::vector<arma::vec> &mus){
  struct Func{
    static inline arma::vec& get_elem(std::vector<arma::vec>::iterator &it){
      return *it;
    }
  };

  return(
    compute_mu_n_Sigma_from_normal_apprx_w_particles
    <densities, std::vector<arma::vec>::iterator, Func, is_forward, use_prior_in_bw_filter>
    (data, t, Q, mus.begin(), mus.size()));
}

/* ------------------------------------------- */

/*
 Output class for smoothers
*/

struct smoother_output {
  struct pair {
    const particle *p;
    double log_weight;

    // added to be able to use std::Vector<>::emplace_back
    pair(particle *p, double log_weight):
      p(p), log_weight(log_weight) {}

    pair():
      p(nullptr),
      log_weight(std::numeric_limits<double>::quiet_NaN()) {}
  };

  struct particle_pairs {
    std::vector<pair> transition_pairs;
    const particle *p;
  };

  std::vector<cloud> forward_clouds;
  std::vector<cloud> backward_clouds;
  std::vector<cloud> smoothed_clouds;
  std::vector<std::vector<particle_pairs>> transition_likelihoods;
};

Rcpp::List get_rcpp_list_from_cloud(
    const smoother_output &sm_output, const PF_data *data = nullptr);

Rcpp::List get_rcpp_list_from_cloud(
    const std::vector<cloud> &clouds, const bool reverse,
    unsigned int state_dim, const PF_data *data = nullptr);

smoother_output get_clouds_from_rcpp_list(const Rcpp::List &rcpp_list);

#undef USE_PRIOR_IN_BW_FILTER_DEFAULT
#undef MAX
#endif
