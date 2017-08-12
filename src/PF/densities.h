#ifndef DENSITIES
#define DENSITIES

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#include "PF_data.h"
#include "particles.h"
#include "dmvnrm.h"

/*
  Each class has the following static functions:
    log_prob_y_given_state:
      Computes log(P(y | alpha_t))
    log_prob_state_given_previous
      Computes log(P(alpha_t | alpha_{t - 1}))
    log_prob_state_given_next
      Computes log(P(alpha_t | alpha_{t + 1}))
    log_prob_state_given_both
      Computes log(P(alpha_t | alpha_{t - 1}, alpha_{t + 1}))
    log_artificial_prior
      Returns the log density of of artificial prior gamma_t(alpha_t)
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

  const arma::mat& get_Q_t(const PF_data &data, arma::uword t){
    if(Q_t_chol_inv.find(t) == Q_t_chol_inv.end())
      Q_t_chol_inv.insert(std::make_pair(t, data.Q.chol_inv / sqrt(t)));

    return Q_t_chol_inv[t];
  }

public:
  static double log_prob_y_given_state(
      const PF_data &data, const particle &p, int t){
    return(log_prob_y_given_state(data, p.state, t));
  }

  static double log_prob_y_given_state(
      const PF_data &data, const arma::vec &coefs, int t){
    const arma::uvec r_set = Rcpp::as<arma::uvec>(data.risk_sets[t - 1]) - 1;
    arma::vec eta =  coefs.t() * data.X.cols(r_set);
    const arma::uvec is_event = data.is_event_in_bin(r_set) == t - 1; /* zero indexed while t is not */

    auto it_eta = eta.begin();
    auto it_is_event = is_event.begin();
    auto n = eta.n_elem;

    double log_like = 0;
    for(arma::uword i = 0; i < n; ++i, ++it_eta, ++it_is_event){
      /* Truncate */
      *it_eta = MIN(MAX(*it_eta, -15), 15);

      log_like += (*it_is_event == 1) ?
        log(1 / (1 + exp(-*it_eta))) : log(1 - 1 / (1 + exp(-*it_eta)));
    }

    if(data.debug > 4){
      data.log(5) << "Computing log(P(y_t|alpha_t)) at time " << t
                  << " with " << eta.n_elem << " observations. "
                  << "The log likelihood is " << log_like
                  << " and the state is:";
      data.log(5) << coefs.t();
    }

    return log_like;
  }

  static double log_prob_state_given_previous(
      const PF_data &data, const particle &p, int t){
    return(dmvnrm_log(p.state, p.parent->state, data.Q.chol_inv));
  }

  static double log_prob_state_given_next(
      const PF_data &data, const particle &p, int t){
    return(dmvnrm_log(p.state, p.parent->state, data.Q.chol_inv));
  }

  double log_artificial_prior(
      const PF_data &data, const particle &p, int t){
    return(dmvnrm_log(p.state, data.a_0 /* note a_o */, get_Q_t(data, t)));
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
