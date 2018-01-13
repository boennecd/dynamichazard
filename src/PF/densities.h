#ifndef DENSITIES
#define DENSITIES

#include "PF_data.h"
#include "particles.h"
#include "dmvnrm.h"
#include "PF_utils.h"
#include "../family.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*
  Each class has the following static functions:
    log_prob_y_given_state:
      Computes log(P(y | alpha_t))
    log_prob_state_given_previous
      Computes log(P(alpha_t | alpha_{t - 1}))
    log_prob_state_given_next
      Computes log(P(alpha_t | alpha_{t + 1}))
    log_artificial_prior
      Returns the log density of of artificial prior gamma_t(alpha_t)
    get_artificial_prior_covar
      Returns the artificial prior covariance matrix at time t
    get_artificial_prior_covar_state_dim
      Returns
        P_t = \left\{
           \begin{matrix}
              R Q R^\top & t = 0 \\
              F P_{t-1} F^\top + R Q R^\top & t > 0
           \end{matrix}
        \right.
    get_artificial_prior_mean
      Returns the artificial prior mean at time t
    get_artificial_prior_mean_state_dim
      Returns m_t = F^t a_0
    d_log_like
      Returns the first deriative of the log likelihood given outcome and the
      state for a given individual
    dd_log_like
      Same as d_log_like for the second derivative
*/

/*
  base densinty class w/ template density functions
*/
template<class outcome_dens>
class pf_dens_base {
  using uword = arma::uword;
  using time_mat_map = std::map<uword /* time */, const arma::mat>;

  /* Maybe checkout this pot regarding thread safty of std::map
   *   https://stackoverflow.com/q/8234633/5861244
   * As of 11/01/2018 all threads will only need the same element t at the same
   * time. Thus, we just need a lock so the first threat adds the element to
   * start with and we should not have any issues from that point of         */
  time_mat_map Q_t_chol_inv;
  time_mat_map P_t;
  std::map<uword /* time */, const arma::vec> m_t;

  const PF_data &data_;

  // TODO: change to pointers to avoid copies
  const arma::mat& get_Q_t_chol_inv(uword t){
#ifdef _OPENMP
#pragma omp critical(get_Q_t_chol_inv_lock)
{
#endif
  if(Q_t_chol_inv.find(t) == Q_t_chol_inv.end()){
    Q_t_chol_inv.insert(std::make_pair(
        t,
        arma::inv(
          arma::trimatu(
            arma::chol(
              get_artificial_prior_covar(t))))));
  }
#ifdef _OPENMP
}
#endif

    return Q_t_chol_inv[t];
  }

public:
  pf_dens_base(const PF_data &data): data_(data) {};

  /* static member functions */
  static double log_prob_y_given_state(
      const PF_data &data, const particle &p, int t){
    return(log_prob_y_given_state(data, p.get_state(), t));
  }

  static double log_prob_y_given_state(
      const PF_data &data, const arma::vec &state,
      int t, arma::uvec &r_set, const bool multithreaded = true){
    const arma::vec coefs = data.lp_map(state).sv;

    double bin_start, bin_stop;
    if(outcome_dens::uses_at_risk_length){
      auto tmp = get_bin_times(data, t);
      bin_start = tmp.start;
      bin_stop = tmp.stop;
    } else
      // avoid wmaybe-uninitialized
      bin_start = bin_stop = std::numeric_limits<double>::quiet_NaN();

    auto jobs = get_work_blocks(r_set.begin(), r_set.end(), data.work_block_size);
    unsigned int n_jobs = jobs.size();

    /* compute log likelihood */
    double log_like = 0;
#ifdef _OPENMP
    /*
      Use lock as critical section will not do if this function is called in
      nested parallel setup. See https://stackoverflow.com/a/20447843
    */
    omp_lock_t *lock;
    if(multithreaded){
      lock = new omp_lock_t;
      omp_init_lock(lock);
    } else
      // avoid wmaybe-uninitialized
      lock = nullptr;
#pragma omp parallel if(multithreaded)
{
#endif
    double my_log_like = 0;
    arma::vec starts;
    arma::vec stops;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned int i = 0; i < n_jobs; ++i){
      auto &job = jobs[i];
      arma::uvec my_r_set(job.start, job.block_size, false /* don't copy */);
      arma::vec eta =  get_linear_product(coefs, data.X, my_r_set);
      if(data.any_fixed_in_M_step)
        eta +=
          get_linear_product(data.fixed_parems, data.fixed_terms, my_r_set);

      const arma::uvec is_event = data.is_event_in_bin(my_r_set) == t - 1; /* zero indexed while t is not */

      if(outcome_dens::uses_at_risk_length){
        starts = data.tstart(my_r_set);
        stops = data.tstop(my_r_set);
      }

      auto it_eta = eta.begin();
      auto it_is_event = is_event.begin();
      auto it_start = starts.begin();
      auto it_stops = stops.begin();
      unsigned int n = eta.n_elem;
      for(uword i = 0; i < n; ++i, ++it_eta, ++it_is_event){
        double at_risk_length = 0;
        if(outcome_dens::uses_at_risk_length){
          at_risk_length = get_at_risk_length(
            *(it_stops++) /* increament here */, bin_stop,
            *(it_start++) /* increament here */, bin_start);
        }

        auto trunc_eta = outcome_dens::truncate_eta(
          *it_is_event, *it_eta, exp(*it_eta), at_risk_length);
        my_log_like += outcome_dens::log_like(
          *it_is_event, trunc_eta, at_risk_length);
      }
    }

#ifdef _OPENMP
    if(multithreaded)
      omp_set_lock(lock);
#endif

    log_like += my_log_like;

#ifdef _OPENMP
    if(multithreaded)
      omp_unset_lock(lock);
#endif
#ifdef _OPENMP
}
    if(multithreaded){
      omp_destroy_lock(lock);
      delete lock;
    }
#endif

    if(data.debug > 4){
      data.log(5) << "Computing log(P(y_t|alpha_t)) at time " << t
                  << " with " << r_set.n_elem << " observations. "
                  << "The log likelihood is " << log_like
                  << " and the state is:" << std::endl
                  << state.t();
    }

    return log_like;
  }

  static double
  log_prob_y_given_state(
    const PF_data &data, const arma::vec &state,
    int t, const bool multithreaded = true){
    arma::uvec r_set = get_risk_set(data, t);

    return log_prob_y_given_state(data, state, t, r_set, multithreaded);
  }

  static double log_prob_state_given_previous(
      const PF_data &data, const arma::vec state,
      arma::vec previous_state, int t){
    return(dmvnrm_log(
        data.err_state_map_inv(state).sv,
        data.err_state_map_inv(
          data.state_trans_map    (previous_state).sv).sv,
        data.Q.chol_inv));
  }

  static double log_prob_state_given_next(
      const PF_data &data, const arma::vec state,
      arma::vec next_state, int t){
    return(dmvnrm_log(
        data.err_state_map_inv(state).sv,
        data.err_state_map_inv(
          data.state_trans_map_inv(    next_state).sv).sv,
        data.Q.chol_inv));
  }

  /* non-static member functions*/
  double log_artificial_prior(const particle &p, int t){
    return(dmvnrm_log(
        data_.err_state_map_inv(p.get_state()).sv,
        get_artificial_prior_mean(t),
        get_Q_t_chol_inv(t)));
  }

  // TODO: change to pointers to avoid copies
  const arma::mat& get_artificial_prior_covar_state_dim(uword t){
#ifdef _OPENMP
#pragma omp critical(get_P_t_lock)
{
#endif
  if(P_t.find(t) == P_t.end()){
    arma::mat new_terms(data_.err_state_map(data_.Q.mat).sv);
    arma::mat out(new_terms);
    for(uword i = 1; i <= t; ++i){
      out = new_terms + data_.state_trans_map(out).sv;

      P_t.insert(std::make_pair(i, out));
    }
  }
#ifdef _OPENMP
}
#endif

    return P_t[t];
  }

  arma::mat get_artificial_prior_covar(uword t){
    return data_.err_state_map_inv(
      get_artificial_prior_covar_state_dim(t)).sv;
  }

  // TODO: change to pointers to avoid copies
  const arma::vec& get_artificial_prior_mean_state_dim(uword t){
#ifdef _OPENMP
#pragma omp critical(get_m_t_lock)
{
#endif
  if(m_t.find(t) == m_t.end()){
    arma::vec out(data_.a_0);
    for(uword i = 1; i <= t; ++i){
      out = data_.state_trans_map(out).sv;

      m_t.insert(std::make_pair(i, out));
    }
  }
#ifdef _OPENMP
}
#endif

    return m_t[t];
  }

  arma::vec get_artificial_prior_mean(uword t){
    return data_.err_state_map_inv(get_artificial_prior_mean_state_dim(t)).sv;
  }
};

/*
  Class for binary outcomes with the logistic link function where state vectors
  follows a first order random walk
*/

class logistic_dens :
  public pf_dens_base<logistic>,
  public logistic {
public:
  using pf_dens_base<logistic>::pf_dens_base;
};

/*
  Class for exponentially distributed arrival times where state vectors follows
  a first order random walk
*/
class exponential_dens :
  public pf_dens_base<exponential>,
  public exponential {
public:
  using pf_dens_base<exponential>::pf_dens_base;
};

#endif
