#ifndef IMPORTANCE_SAMPLERS
#define IMPORTANCE_SAMPLERS

#define MAX(a,b) (((a)>(b))?(a):(b))

#include "PF_data.h"
#include "particles.h"
#include "../sample_funcs.h"
#include "dmvnrm.h"
#include "PF_utils.h"
#include "../BLAS_and_LAPACK/arma_utils.h"

/*
  Each "importance_sampler" has two static functions:
    sample:
      Returns a new cloud of particles with sampled states. The states are
      sampled acording to the specific importance density. The function also
      sets log_importance_dens on each particle
    sample_smooth:
      Returns sample given past and next state
    log_importance_dens_smooth:
      Returns the log importance density for state smoothing
    sample_first_state_n_set_weights:
      Returns a particle cloud for time zero or d + 1 with weights set
*/

/* base class importance samplers */

template<typename densities, bool is_forward>
class importance_dens_base {
public:
  static cloud sample_first_state_n_set_weights(const PF_data &data){
    cloud ans;

    if(is_forward){
      ans.reserve(data.N_fw_n_bw);
      double log_weight = log(1. / data.N_fw_n_bw);
      for(arma::uword i = 0; i < data.N_fw_n_bw; ++i) {
        // Do not sample but set to a_0 with equal weight on each
        ans.New_particle(data.a_0, nullptr);
        ans[i].log_weight = log_weight;
      }

    } else {
      ans.reserve(data.N_first);
      const arma::mat Q_d_chol = data.Q.chol * data.d + 1;
      const arma::mat Q_chol_inv = arma::inv(arma::trimatu(Q_d_chol));
      double max_weight =  -std::numeric_limits<double>::max();
      for(arma::uword i = 0; i < data.N_first; ++i){
        arma::vec new_state = mvrnorm(data.a_0, Q_d_chol);
        ans.New_particle(new_state, nullptr);
        ans[i].log_weight = dmvnrm_log(new_state, data.a_0, Q_chol_inv);

        max_weight = MAX(ans[i].log_weight, max_weight);
      }

      normalize_log_weights<false, true>(ans, max_weight);
    }

    return(ans);
  }
};

/*
 Importance sampler with importance density that does not depend on the
 outcome for the first order random walk. That is:
  q(alpah_t | alpha_{t - 1}^{j}, y_t) =
    N(alpha | alpha_{t - j}^{j}, Q)
  tilde{q}(alpah_t | alpha_{t + 1}^{k}, y_t) =
    N(alpha | alpha_{t + 1}^{j}, Q)
  tilde{q}(alpah_t | alpha_{t + 1}^{j}, alpha_{t + 1}^{k}, y_t) =
    N(alpha | (1/2)Q(Q^-1 alpha_{t + 1}^{j} + alpha_{t + 1}^{k}), (1/2)Q) =
    N(alpha | (1/2)alpha_{t + 1}^{j} + (1/2)alpha_{t + 1}^{k}), (1/2)Q)

 See:
  Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class importance_dens_no_y_dependence :
  public importance_dens_base<densities, is_forward>{
  static double log_importance_dens(const PF_data &data, const particle &p, int t){
    /* independent of is_forward for first order random walk */
    return dmvnrm_log(p.state, p.parent->state, data.Q_proposal.chol_inv);
  }

  static double log_importance_dens_smooth(
      const PF_data &data, const particle &p, int t){
    arma::vec mean = p.parent->state + p.child->state;
    mean *= .5;

    return dmvnrm_log(p.state, mean, data.Q_proposal_smooth.chol_inv);
  }

public:
  static cloud sample(const PF_data &data, cloud &cl, const arma::uvec &resample_idx, const unsigned int t){
    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by";
      data.log(3) << data.Q_proposal.chol;
    }

    auto it = resample_idx.begin();
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      arma::vec new_state = mvrnorm(cl[*it].state, data.Q_proposal.chol);
      ans.New_particle(new_state, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = log_importance_dens(data, p, t);
    }

    return ans;
  }

  static cloud sample_smooth(
    const PF_data &data,
    cloud &fw_cloud, const arma::uvec &fw_idx,
    cloud &bw_cloud, const arma::uvec &bw_idx, const unsigned int t){
    cloud ans;
    ans.reserve(data.N_smooth);

    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Q) given by";
      data.log(3) << data.Q_proposal_smooth.chol;
    }

    auto it_fw = fw_idx.begin();
    auto it_bw = bw_idx.begin();
    for(arma::uword i = 0; i < data.N_smooth; ++i, ++it_fw, ++it_bw){
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mean = fw_p.state + bw_p.state;
      mean *= .5;

      arma::vec new_state = mvrnorm(
        mean, data.Q_proposal_smooth.chol);
      ans.New_particle(new_state, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = log_importance_dens_smooth(data, p, t);
    }

    return ans;
  }
};

/*
  Sampler which makes a normal approximation for the observed outcome. See the
  second example on page 462-463 of:
    Fearnhead, P., Wyncoll, D., & Tawn, J. (2010). A sequential smoothing algorithm with linear computational cost. Biometrika, 97(2), 447-464.
*/

template<typename densities, bool is_forward>
class importance_dens_normal_approx  :
  public importance_dens_base<densities, is_forward>{
  struct intermediate_output {
    arma::vec mu;
    arma::mat Sigma_inv_chol;
    arma::mat Sigma_chol;
    arma::mat sigma_chol_inv;
  };

  static intermediate_output compute_mu_n_Sigma(
      const PF_data &data, const unsigned int t, const PF_data::covarmat &Q, const arma::vec &alpha_bar){
    if(data.debug > 2){
      data.log(3) << "Computing normal approximation with mean vector:";
      data.log(3) << alpha_bar.t();
      data.log(3) << "and covaraince matrix:";
      data.log(3) << Q.mat;
    }

    intermediate_output ans;

    /* Compute the terms that does not depend on the outcome */
    /* Sigma^-1 = (Q + \tilde{Q})^{-1} */
    arma::mat Sigma_inv = Q.inv;
    ans.mu = arma::vec(alpha_bar.n_elem, arma::fill::zeros);
    arma::vec &mu = ans.mu;

    /* Add the terms that does depend on the outcome */
    const arma::uvec r_set = Rcpp::as<arma::uvec>(data.risk_sets[t - 1]) - 1;
    arma::vec eta =  alpha_bar.t() * data.X.cols(r_set);
    const arma::uvec is_event = data.is_event_in_bin(r_set) == t - 1; /* zero indexed while t is not */

    auto it_eta = eta.begin();
    auto it_is_event = is_event.begin();
    auto it_r = r_set.begin();
    arma::uword n_elem = eta.n_elem;
    /*
      Update with:
      Signa = ... + X^T (-G) X
      mu = X^T (-G) X \bar{alpha} + X^T (-g)
    */
    for(arma::uword i = 0; i < n_elem; ++i, ++it_eta, ++it_is_event, ++it_r){
      double g = densities::log_p_prime(*it_is_event, *it_eta, t);
      double neg_G = - densities::log_p_2prime(*it_is_event, *it_eta, t);

      sym_mat_rank_one_update(neg_G, data.X.col(*it_r), Sigma_inv);

      mu += data.X.col(*it_r) * ((*it_eta * neg_G) + g);
    }

    /* copy to lower */
    Sigma_inv = arma::symmatu(Sigma_inv);

    /* Compute needed factorizations */
    ans.Sigma_inv_chol = arma::chol(Sigma_inv);
    ans.Sigma_chol = arma::chol(arma::inv(Sigma_inv)); // TODO: do something smarter
    ans.sigma_chol_inv = arma::inv(arma::trimatu(ans.Sigma_chol));

    return ans;
  }

  inline static void debug_msg_before_sampling(const PF_data &data, const intermediate_output &inter_output){
    if(data.debug > 2){
      data.log(3) << "Sampling new cloud from normal distribution with chol(Sigma) given by";
      data.log(3) << inter_output.Sigma_chol;
      data.log(3) << "The mean before accounting for the parent (and child) particle is:";
      data.log(3) << solve_w_precomputed_chol(inter_output.Sigma_inv_chol, inter_output.mu).t();
    }
  }

  inline static void debug_msg_while_sampling(
      const PF_data &data, const particle &p, const arma::vec &mu){
    if(data.debug > 4){
      data.log(5) << "Sampled particle:";
      data.log(5) << p.state.t();
      data.log(5) << "from normal distribution with mean:";
      data.log(5) << mu.t();
      data.log(5) << "The parent had state:";
      data.log(5) << p.parent->state.t();

      if(p.child){
        data.log(5) << "and the child had state";
        data.log(5) << p.child->state.t();
      }
    }
  }

public:
  static cloud sample(const PF_data &data, cloud &cl, const arma::uvec &resample_idx, const unsigned int t){
    /* Find weighted mean estimate */
    arma::uword p = cl[0].state.n_elem;
    arma::vec alpha_bar(p, arma::fill::zeros);
    for(auto it = cl.begin(); it != cl.end(); ++it){
      alpha_bar += (it->state * exp(it->log_weight));
    }

    /* compute parts of the terms for the mean and covariance */
    auto &Q = data.Q_proposal;
    auto inter_output = compute_mu_n_Sigma(data, t, Q, alpha_bar);

    /* Sample */
    debug_msg_before_sampling(data, inter_output);

    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    auto it = resample_idx.begin();
    for(arma::uword i = 0; i < data.N_fw_n_bw; ++i, ++it){
      arma::vec mu_j = solve_w_precomputed_chol(Q.chol, cl[*it].state) + inter_output.mu;
      mu_j = solve_w_precomputed_chol(inter_output.Sigma_inv_chol, mu_j);

      arma::vec new_state = mvrnorm(mu_j, inter_output.Sigma_chol);
      ans.New_particle(new_state, &cl[*it]);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.state, mu_j, inter_output.sigma_chol_inv);

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }

  static cloud sample_smooth(
      const PF_data &data,
      cloud &fw_cloud, const arma::uvec &fw_idx,
      cloud &bw_cloud, const arma::uvec &bw_idx, const unsigned int t){
    /* Find weighted mean estimate */
    arma::uword p = fw_cloud[0].state.n_elem;
    arma::vec alpha_bar(p, arma::fill::zeros);
    for(auto it = fw_cloud.begin(); it != fw_cloud.end(); ++it){
      alpha_bar += (it->state * (.5 * exp(it->log_weight)));
    }
    for(auto it = bw_cloud.begin(); it != bw_cloud.end(); ++it){
      alpha_bar += (it->state * (.5 * exp(it->log_weight)));
    }

    /* compute parts of the terms for the mean and covariance */
    auto &Q = data.Q_proposal_smooth;
    auto inter_output = compute_mu_n_Sigma(data, t, Q, alpha_bar);

    /* Sample */
    debug_msg_before_sampling(data, inter_output);

    cloud ans;
    ans.reserve(data.N_fw_n_bw);

    auto it_fw = fw_idx.begin();
    auto it_bw = bw_idx.begin();
    for(arma::uword i = 0; i < data.N_smooth; ++i, ++it_fw, ++it_bw){
      const particle &fw_p = fw_cloud[*it_fw];
      const particle &bw_p = bw_cloud[*it_bw];

      arma::vec mu_j = fw_p.state + bw_p.state;
      mu_j *= .5;
      mu_j = solve_w_precomputed_chol(Q.chol, mu_j) + inter_output.mu;
      mu_j = solve_w_precomputed_chol(inter_output.Sigma_inv_chol, mu_j);

      arma::vec new_state = mvrnorm(mu_j, inter_output.Sigma_chol);
      ans.New_particle(new_state, &fw_p, &bw_p);

      particle &p = ans[i];
      p.log_importance_dens = dmvnrm_log(p.state, mu_j, inter_output.sigma_chol_inv);

      debug_msg_while_sampling(data, p, mu_j);
    }

    return(ans);
  }
};



#undef MAX
#endif
