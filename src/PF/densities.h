#ifndef DENSITIES
#define DENSITIES

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#include "PF_data.h"
#include "particles.h"

/*
  Each class has the following static functions:
    log_prob_y_given_state:
      Computes log(P(y | alpha_t))
    log_prob_state_given_previous
      Computes log(P(alpha_t | alpha_{t - 1}))
    log_prob_state_given_next
      Computes log(P(alpha_t | alpha_{t + 1}))
*/

class binary {
public:
  static double log_prob_y_given_state(
      const PF_data &data, particle &p, int t){
    const arma::uvec r_set = Rcpp::as<arma::uvec>(data.risk_sets[t - 1]) - 1;
    arma::vec eta = p.state * data.X.cols(r_set);
    const arma::uvec is_event = data.is_event_in_bin(r_set) == t;

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

    return log_like;
  }
};

#endif
