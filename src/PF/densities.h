#ifndef DENSITIES
#define DENSITIES

#include "dmvnrm.h"
#include "PF_data.h"
#include "get_work_blocks.h"
#include "../family.h"
#include "particles.h"

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
    d_log_like
      Returns the first deriative of the log likelihood given outcome and the
      state for a given individual
    dd_log_like
      Same as d_log_like for the second derivative
*/

class pf_base_dens : public virtual family_base {
  using uword = arma::uword;

protected:
  const PF_data &data;

public:
  pf_base_dens(const PF_data &data): data(data) {};
  virtual ~pf_base_dens() = default;

  double log_prob_state_given_previous(
      const arma::vec state, arma::vec previous_state, int t){
    return(dmvnrm_log(
        data.err_state_inv->map(state).sv,
        data.err_state_inv->map(
            data.state_trans->map(previous_state).sv).sv,
            data.Q.chol_inv()));
  }

  double log_artificial_prior(const particle &p, int t){
    return(dmvnrm_log(
        p.get_state(),
        data.uncond_mean_state(t),
        data.uncond_covar_state(t).chol_inv()));
  }

  double log_prob_y_given_state(
      const arma::vec &state, int t, arma::uvec &r_set,
      const bool multithreaded = true) const {
    const arma::vec coefs = data.err_state_inv->map(state).sv;
    bool uses_at_risk_len = uses_at_risk_length();

    double bin_start, bin_stop;
    if(uses_at_risk_len){
      auto tmp = get_bin_times(data, t);
      bin_start = tmp.start;
      bin_stop = tmp.stop;
    } else
      // avoid wmaybe-uninitialized
      bin_start = bin_stop = std::numeric_limits<double>::quiet_NaN();

    auto jobs = get_work_blocks(
      r_set.begin(), r_set.end(), data.work_block_size);
    unsigned int n_jobs = jobs.size();

    /* compute log likelihood */
    double result = 0;
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
    double my_result = 0;
    arma::vec starts;
    arma::vec stops;

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for(unsigned int i = 0; i < n_jobs; ++i){
      auto &job = jobs[i];
      arma::uvec my_r_set(job.start, job.block_size, false /* don't copy */);
      arma::vec eta =  get_linear_product(coefs, data.X, my_r_set);
      eta += data.fixed_effects(my_r_set);

      const arma::uvec is_event = data.is_event_in_bin(my_r_set) == t - 1; /* zero indexed while t is not */

          if(uses_at_risk_len){
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
            if(uses_at_risk_len){
              at_risk_length = get_at_risk_length(
                *(it_stops++) /* increament here */, bin_stop,
                *(it_start++) /* increament here */, bin_start);
            }

            auto trunc_eta = truncate_eta(
              *it_is_event, *it_eta, exp(*it_eta), at_risk_length);
            my_result += log_like(
              *it_is_event, trunc_eta, at_risk_length);
          }
    }

#ifdef _OPENMP
  if(multithreaded)
    omp_set_lock(lock);
#endif

    result += my_result;

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
                  << "The log likelihood is " << result
                  << " and the state is:" << std::endl
                  << state.t();
    }

    return result;
  }

  double
    log_prob_y_given_state(
      const arma::vec &state, int t, const bool multithreaded = true){
      arma::uvec r_set = get_risk_set(data, t);

      return log_prob_y_given_state(state, t, r_set, multithreaded);
    }

  double log_prob_y_given_state(const particle &p, int t){
    return(log_prob_y_given_state(p.get_state(), t));
  }
};

/* Class for binary outcomes with logistic link */

class logistic_dens :
  public virtual pf_base_dens,
  public virtual logistic {
public:
  logistic_dens(const PF_data &data): pf_base_dens(data) {};
};

/* Class for piece-wise constant exponentially distributed arrival times */
class exponential_dens :
  public virtual pf_base_dens,
  public virtual exponential {
public:
  exponential_dens(const PF_data &data): pf_base_dens(data) {};
};

#endif
