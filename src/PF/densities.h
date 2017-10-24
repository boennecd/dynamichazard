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
    log_artificial_prior
      Returns the log density of of artificial prior gamma_t(alpha_t)
    get_artificial_prior_covar
      Returns the artificial prior covariance matrix at time t
    d_log_like
      Returns the first deriative of the log likelihood given outcome and the
      state for a given individual
    dd_log_like
      Same as d_log_like for the second derivative
*/

/*
  Class for binary outcomes with the logistic link function where state vectors
  follows a first order random walk
*/
template<class outcome_dens>
class first_order_random_walk_base {
  using uword = arma::uword;
  using chol_map = std::map<uword /* time */, const arma::mat /* inv chol cov matrix */>;
  chol_map Q_t_chol_inv;

  const arma::mat& get_Q_t_chol_inv(const PF_data &data, arma::uword t){
#ifdef _OPENMP
#pragma omp critical(get_Q_t_chol_inv_lock)
{
#endif
    if(Q_t_chol_inv.find(t) == Q_t_chol_inv.end()){
      Q_t_chol_inv.insert(std::make_pair(
          t,
          arma::inv(arma::trimatu(arma::chol(get_artificial_prior_covar(
            data, t))))));
    }
#ifdef _OPENMP
}
#endif

    return Q_t_chol_inv[t];
  }

public:
  static double log_prob_y_given_state(
      const PF_data &data, const particle &p, int t){
    return(log_prob_y_given_state(data, p.get_state(), t));
  }

  static double log_prob_y_given_state
  (const PF_data &data, const arma::vec &coefs,
   int t, arma::uvec &r_set, const bool multithreaded = true){
    double bin_start, bin_stop;
    if(outcome_dens::uses_at_risk_length){
      auto tmp = get_bin_times(data, t);
      bin_start = tmp.start;
      bin_stop = tmp.stop;
    }

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
    }
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
      for(arma::uword i = 0; i < n; ++i, ++it_eta, ++it_is_event){
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
                  << coefs.t();
    }

    return log_like;
  }

  static double log_prob_y_given_state
  (const PF_data &data, const arma::vec &coefs,
   int t, const bool multithreaded = true){
    arma::uvec r_set = get_risk_set(data, t);

    return log_prob_y_given_state(data, coefs, t, r_set, multithreaded);
  }

  static double log_prob_state_given_previous(
      const PF_data &data, const arma::vec state, const arma::vec previous_state, int t){
    return(dmvnrm_log(state, previous_state, data.Q.chol_inv));
  }

  double log_artificial_prior(
      const PF_data &data, const particle &p, int t){
    return(dmvnrm_log(p.get_state(), data.a_0 /* note a_o */, get_Q_t_chol_inv(data, t)));
  }

  static arma::mat get_artificial_prior_covar(
      const PF_data &data, int t){
    return data.Q.mat * (t + 1);
  }
};

/*
  Class for binary outcomes with the logistic link function where state vectors
  follows a first order random walk
*/

class logistic_dens :
  public first_order_random_walk_base<logistic>,
  public logistic {};

/*
  Class for exponentially distributed arrival times where state vectors follows
  a first order random walk
*/
class exponential_dens :
  public first_order_random_walk_base<exponential>,
  public exponential {};

#endif
