#ifndef DENSITIES
#define DENSITIES

#include "PF_data.h"
#include "particles.h"
#include "dmvnrm.h"
#include "PF_utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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
    log_p_prime
      Returns the first deriative of the log likelihood given outcome and the
      state for a given individual
    log_p_2prime
      Same as log_p_prime for the second derivative
*/

/*
  Class for binary outcomes with the logistic link function where state vectors
  follows a first order random walk
*/
class binary {
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
    return(log_prob_y_given_state(data, p.state, t));
  }

  static double log_prob_y_given_state
  (const PF_data &data, const arma::vec &coefs,
   int t, arma::uvec &r_set, const bool multithreaded = true){
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
#pragma omp parallel for schedule(static) if(multithreaded)
#endif
    for(unsigned int i = 0; i < n_jobs; ++i){
      auto &job = jobs[i];
      arma::uvec my_r_set(job.start, job.block_size, false /* don't copy */);
      arma::vec eta =  coefs.t() * data.X.cols(my_r_set);
      const arma::uvec is_event = data.is_event_in_bin(my_r_set) == t - 1; /* zero indexed while t is not */

      auto it_eta = eta.begin();
      auto it_is_event = is_event.begin();
      auto n = eta.n_elem;
      double my_log_like = 0;
      for(arma::uword i = 0; i < n; ++i, ++it_eta, ++it_is_event){
        /* Truncate */
        *it_eta = MIN(MAX(*it_eta, -15), 15);

        my_log_like += (*it_is_event == 1) ?
          log(1 / (1 + exp(-*it_eta))) : log(1 - 1 / (1 + exp(-*it_eta)));
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
    }
#ifdef _OPENMP
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
    arma::uvec r_set = data.get_risk_set(t);

    return log_prob_y_given_state(data, coefs, t, r_set, multithreaded);
  }

  static double log_prob_state_given_previous(
      const PF_data &data, const arma::vec state, const arma::vec previous_state, int t){
    return(dmvnrm_log(state, previous_state, data.Q.chol_inv));
  }

  double log_artificial_prior(
      const PF_data &data, const particle &p, int t){
    return(dmvnrm_log(p.state, data.a_0 /* note a_o */, get_Q_t_chol_inv(data, t)));
  }

  static arma::mat get_artificial_prior_covar(
      const PF_data &data, int t){
    return data.Q.mat * (t + 1);
  }

  static double log_p_prime(double y, double eta, int t){
    if(eta < 0){
      double exp_eta = exp(eta);
      return (exp_eta * (y - 1) + y) / (exp_eta + 1);
    }

    double exp_neg_eta = exp(-eta);
    return ((y - 1) + y * exp_neg_eta) / (1 + exp_neg_eta);
  }

  static double log_p_2prime(double y, double eta, int t){
    if(eta < 0){
      double exp_eta = exp(eta);
      return - exp_eta / ((1 + exp_eta) * (1 + exp_eta));
    }

    double exp_neg_eta = exp(-eta);
    return - exp_neg_eta / ((1 + exp_neg_eta) * (1 + exp_neg_eta));
  }
};

#undef MIN
#undef MAX
#endif
