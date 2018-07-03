#include "PF_utils.h"

/* Perform normal approximation through second order Taylor approximation as
 * in
 *     Pitt, M. K., & Shephard, N. (1999). Filtering via simulation: Auxiliary
 *    particle filters. Journal of the American statistical association,
 *    94(446), 590-599.
 *
 * page 594.
 */

input_for_normal_apprx taylor_normal_approx(
    pf_base_dens &dens_calc, const PF_data &data,
    const unsigned int t, const covarmat &Q, const arma::vec &alpha_bar,
    arma::uvec &r_set, const unsigned int debug_lvl,
    const bool multithread){
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
  } // end omp parallel
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

input_for_normal_apprx taylor_normal_approx(
    pf_base_dens &dens_calc, const PF_data &data, const unsigned int t,
    const covarmat &Q, const arma::vec &alpha_bar,
    const unsigned int debug_lvl, const bool multithread){
  /*
   Had similar issues as posted here:
   http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-June/005968.html

   Thus, I made this overload
   */

  arma::uvec r_set = get_risk_set(data, t);

  return(
    taylor_normal_approx
    (dens_calc, data, t, Q, alpha_bar, r_set, debug_lvl, multithread));
}

/*----------------------------------------------------*/

input_for_normal_apprx_w_cloud_mean
  taylor_normal_approx_w_cloud_mean(
    pf_base_dens &dens_calc, const PF_data &data,
    const unsigned int t, const covarmat &Q, const arma::vec &alpha_bar,
    cloud &cl /* set mu_js when cloud is passed to */, const bool is_forward){
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
      taylor_normal_approx
      (dens_calc, data, t, *Q_use, alpha_bar, 2, true);

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

/*----------------------------------------------------------*/

template<typename mu_iterator, typename Func>
input_for_normal_apprx_w_particle_mean
  taylor_normal_approx_w_particles(
   pf_base_dens &dens_calc, const PF_data &data,
   const unsigned int t,
   const covarmat &Q, mu_iterator begin, const unsigned int size,
   const bool is_forward){
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
   auto inter = taylor_normal_approx
     (dens_calc, data, t, *Q_use, this_state, r_set, 5, false);

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

input_for_normal_apprx_w_particle_mean
  taylor_normal_approx_w_particles(
   pf_base_dens &dens_calc, const PF_data &data,
   const unsigned int t, const covarmat &Q, cloud &cl,
   const bool is_forward){
  struct Func{
    static inline const arma::vec get_elem(cloud::iterator &it){
     return it->get_state();
    }
  };

  return(
    taylor_normal_approx_w_particles
      <cloud::iterator, Func>
      (dens_calc, data, t, Q, cl.begin(), cl.size(), is_forward));
 }

input_for_normal_apprx_w_particle_mean
 taylor_normal_approx_w_particles(
   pf_base_dens &dens_calc, const PF_data &data,
   const unsigned int t, const covarmat &Q, std::vector<arma::vec> &states,
   const bool is_forward){
   struct Func{
     static inline arma::vec& get_elem(std::vector<arma::vec>::iterator &it){
       return *it;
     }
   };

   return(
     taylor_normal_approx_w_particles
     <std::vector<arma::vec>::iterator, Func>
     (dens_calc, data, t, Q, states.begin(), states.size(), is_forward));
 }
