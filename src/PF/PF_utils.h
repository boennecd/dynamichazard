#ifndef PF_UTILS
#define PF_UTILS

#include <tuple>
#include "particles.h"
#include "densities.h"
#include "../arma_BLAS_LAPACK.h"
#include "../utils.h"

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

/* The function below takes in a state \bar{\alpha} and returns
 * mu = R^\top L^\top X^\top ((-G(\bar{\alpha})) X L \bar{\alpha}
 *       + g(\bar{\alpha}))
 * Sigma = (R^\top L^\top X^\top (-G(\bar{\alpha})) X L R + Q^-1)^-1 */

struct input_for_normal_apprx {
  /* mean of error term              */
  arma::vec mu;
  /* covariance matrix of error term */
  arma::mat Sigma_inv;
  arma::mat Sigma_inv_chol;
  arma::mat Sigma_chol;
  arma::mat sigma_chol_inv;
};

template<unsigned int debug_lvl, bool multithread>
static input_for_normal_apprx compute_mu_n_Sigma_from_normal_apprx(
    pf_base_dens &dens_calc,
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
    <debug_lvl, multithread>
    (dens_calc, data, t, Q, alpha_bar, r_set));
}

template<unsigned int debug_lvl, bool multithread>
static input_for_normal_apprx compute_mu_n_Sigma_from_normal_apprx(
    pf_base_dens &dens_calc,
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
  const bool uses_at_risk_length = dens_calc.uses_at_risk_length();
  if(uses_at_risk_length){
    auto tmp = get_bin_times(data, t);
    bin_start = tmp.start;
    bin_stop = tmp.stop;

  } else
    // avoid wmaybe-uninitialized
    bin_start = bin_stop = std::numeric_limits<double>::quiet_NaN();

  /* Compute the terms that does not depend on the outcome */
  /* Sigma^-1 = (Q + \tilde{Q})^{-1} */
  arma::uword p = data.err_dim;
  arma::vec coefs = data.err_state_inv->map(alpha_bar).sv;
  arma::mat Sigma_inv = Q.inv;

  ans.mu = arma::vec(p, arma::fill::zeros);
  arma::vec &mu = ans.mu;

  /* Add the terms that does depend on the outcome */
  auto jobs =
    get_work_blocks(r_set.begin(), r_set.end(), data.work_block_size);
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
  } else
    lock = nullptr;
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

    arma::vec eta =  get_linear_product(coefs, data.X, my_r_set);
    arma::vec offsets = data.fixed_effects(my_r_set);
    eta += offsets;
    const arma::uvec is_event = data.is_event_in_bin(my_r_set) == t - 1; /* zero indexed while t is not */

    if(uses_at_risk_length){
      starts = data.tstart(my_r_set);
      stops = data.tstop(my_r_set);
    }

    auto it_eta = eta.begin();
    auto it_off = offsets.begin();
    auto it_is_event = is_event.begin();
    auto it_r = my_r_set.begin();
    auto it_start = starts.begin();
    auto it_stops = stops.begin();
    arma::uword n_elem = eta.n_elem;
    /*
      Update with:
        Sigma = ... + R^T L^T X^T (-G) X L R
        mu    = R^T L^T X^T (-G) X L \bar{alpha} + R^T L^T X^T (-g)
    */
    for(arma::uword i = 0; i < n_elem;
        ++i, ++it_eta, ++it_is_event, ++it_r, ++it_off){
      double at_risk_length = 0;
      if(uses_at_risk_length){
        at_risk_length = get_at_risk_length(
          *(it_stops++) /* increament here */, bin_stop,
          *(it_start++) /* increament here */, bin_start);

      }

      auto trunc_eta = dens_calc.truncate_eta(
        *it_is_event, *it_eta, exp(*it_eta), at_risk_length);
      double g = dens_calc.d_log_like(
        *it_is_event, trunc_eta, at_risk_length);
      double neg_G = - dens_calc.dd_log_like(
        *it_is_event, trunc_eta, at_risk_length);

      arma::vec x_err_space = data.X.col(*it_r);
      sym_mat_rank_one_update(neg_G, x_err_space, my_Sigma_inv);
      my_mu += x_err_space * (((*it_eta - *it_off) * neg_G) + g);
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
  std::swap(ans.Sigma_inv, Sigma_inv);

  return ans;
}

/* ------------------------------------------- */

struct input_for_normal_apprx_w_cloud_mean : public input_for_normal_apprx {
  /* the conditional means for each of the particles given their parents */
  std::vector<arma::vec> mu_js;

  input_for_normal_apprx_w_cloud_mean(input_for_normal_apprx &&other):
    input_for_normal_apprx(other) {}
};

template<bool is_forward>
static input_for_normal_apprx_w_cloud_mean
  compute_mu_n_Sigma_from_normal_apprx_w_cloud_mean(
    pf_base_dens &dens_calc, const PF_data &data,
    const unsigned int t, const covarmat &Q, const arma::vec &alpha_bar,
    cloud &cl /* set mu_js when cloud is passed to */){
    const covarmat *Q_use;
    const arma::mat tmp;
    const arma::vec *mu_term;
    if(!is_forward){
      Q_use = new covarmat(arma::inv(
        data.state_trans_err_inv->map(Q.inv).sv + data.uncond_covar(t)));
      mu_term = &data.uncond_mean(t);

    } else {
      Q_use = &Q;
      // avoid wmaybe-uninitialized
      mu_term = nullptr;

    }

    input_for_normal_apprx_w_cloud_mean ans =
      compute_mu_n_Sigma_from_normal_apprx
      <2, true>
      (dens_calc, data, t, *Q_use, alpha_bar);

    auto n_elem = cl.size();
    ans.mu_js = std::vector<arma::vec>(n_elem);
#ifdef _OPENMP
#pragma omp  parallel for schedule(static)
#endif
    for(unsigned int i = 0; i < n_elem; ++i){
      arma::vec mu_j;
      if(is_forward){
        mu_j = data.state_trans->map(cl[i].get_state()).sv;
        mu_j = data.err_state_inv->map(mu_j).sv;
        mu_j = solve_w_precomputed_chol(Q.chol, mu_j) + ans.mu;

      } else {
        mu_j = data.err_state_inv->map(cl[i].get_state()).sv;
        mu_j = solve_w_precomputed_chol(Q.chol, mu_j);
        mu_j = data.state_trans_inv->map(data.err_state->map(mu_j).sv).sv;
        mu_j = data.err_state_inv->map(mu_j).sv      + ans.mu + *mu_term;

      }

      mu_j = solve_w_precomputed_chol(ans.Sigma_inv_chol, mu_j);

      ans.mu_js[i] = std::move(mu_j);
    }

    if(!is_forward){
      delete Q_use;

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

template<typename mu_iterator,
         typename Func, bool is_forward>
static input_for_normal_apprx_w_particle_mean
compute_mu_n_Sigma_from_normal_apprx_w_particles(
  pf_base_dens &dens_calc, const PF_data &data,
  const unsigned int t,
  const covarmat &Q, mu_iterator begin, const unsigned int size){
  const covarmat *Q_use;
  arma::mat Q_art_chol;
  const arma::vec *mu_term;
  if(!is_forward){
    // Add the covariance matrix of the artificial prior
    Q_use = new covarmat(arma::inv(Q.inv + data.uncond_covar(t)));
    mu_term = &data.uncond_mean(t);

  } else {
    Q_use = &Q;
    // avoid wmaybe-uninitialized
    mu_term = nullptr;

  }

  input_for_normal_apprx_w_particle_mean ans(size);

  mu_iterator b = begin;
  arma::uvec r_set = get_risk_set(data, t);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(unsigned int i = 0; i < size; ++i){
    mu_iterator iter = b + i;
    const arma::vec &this_state = Func::get_elem(iter);
    auto inter = compute_mu_n_Sigma_from_normal_apprx
      <5, false>(dens_calc, data, t, *Q_use, this_state, r_set);

    arma::vec mu;
    if(is_forward){
      mu = data.state_trans->map(this_state).sv;
      mu = solve_w_precomputed_chol(Q.chol, mu) + inter.mu;

    } else {
      mu = data.err_state_inv->map(this_state).sv;
      mu = solve_w_precomputed_chol(Q.chol, mu);
      mu = data.state_trans_inv->map(data.err_state->map(mu).sv).sv;
      mu = data.err_state_inv->map(mu).sv      + inter.mu + *mu_term;

    }

    mu = solve_w_precomputed_chol(inter.Sigma_inv_chol, mu);

    std::swap(ans[i].mu, mu);
    std::swap(ans[i].sigma_chol_inv, inter.sigma_chol_inv);
    std::swap(ans[i].Sigma_chol, inter.Sigma_chol);
  }

  if(!is_forward){
    delete Q_use;

  }

  return ans;
}

template<bool is_forward>
static input_for_normal_apprx_w_particle_mean
compute_mu_n_Sigma_from_normal_apprx_w_particles(
  pf_base_dens &dens_calc, const PF_data &data,
  const unsigned int t, const covarmat &Q, cloud &cl){
  struct Func{
    static inline const arma::vec get_elem(cloud::iterator &it){
      return it->get_state();
    }
  };

  return(
    compute_mu_n_Sigma_from_normal_apprx_w_particles
    <cloud::iterator, Func, is_forward>
    (dens_calc, data, t, Q, cl.begin(), cl.size()));
}

template<bool is_forward>
static input_for_normal_apprx_w_particle_mean
compute_mu_n_Sigma_from_normal_apprx_w_particles(
  pf_base_dens &dens_calc, const PF_data &data,
  const unsigned int t, const covarmat &Q, std::vector<arma::vec> &states){
  struct Func{
    static inline arma::vec& get_elem(std::vector<arma::vec>::iterator &it){
      return *it;
    }
  };

  return(
    compute_mu_n_Sigma_from_normal_apprx_w_particles
    <std::vector<arma::vec>::iterator, Func, is_forward>
    (dens_calc, data, t, Q, states.begin(), states.size()));
}

/* ------------------------------------------- */

/*
 Output class for smoothers
*/

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

struct PF_summary_stats_RW {
  /* E[X_t] */
  std::vector<arma::vec> E_xs;
  /* E[(X_t - FX_{t-1})(X_t - FX_{t-1})^T] */
  std::vector<arma::mat> E_x_less_x_less_one_outers;
};

PF_summary_stats_RW
  compute_summary_stats_first_o_RW
  (const smoother_output&, const arma::vec&, const arma::mat&,
   const arma::mat&);

#undef USE_PRIOR_IN_BW_FILTER_DEFAULT
#undef MAX
#endif
