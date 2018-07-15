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
    const unsigned int t, const arma::mat &Q_inv, const arma::vec &alpha_bar,
    arma::uvec &r_set, const unsigned int debug_lvl,
    const bool multithread, const bool is_forward){
  if(data.debug > debug_lvl){
    data.log(debug_lvl) << "Computing normal approximation with mean vector:" << std::endl
                        << alpha_bar.t()
                        << "and inv(covariance):"  << std::endl
                        << Q_inv;
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

  /* initalize output */
  const arma::uword r = data.err_dim;
  arma::vec coefs = data.err_state_inv->map(alpha_bar).sv;
  arma::mat Sigma_inv(r, r, arma::fill::zeros);
  ans.mu = arma::vec(r, arma::fill::zeros);
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
  arma::mat my_Sigma_inv(r, r, arma::fill::zeros);
  arma::vec my_mu(r, arma::fill::zeros);

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
    /* Update with inverse of sigma and mu */
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
  if(!is_forward)
    Sigma_inv = data.err_state->map(Sigma_inv).sv;
  Sigma_inv += Q_inv;

  /* Compute needed factorizations */
  ans.Sigma_inv_chol = arma::chol(Sigma_inv);
  ans.Sigma = arma::inv(Sigma_inv);
  ans.Sigma_chol = arma::chol(ans.Sigma);
  ans.sigma_chol_inv = arma::inv(arma::trimatu(ans.Sigma_chol));
  std::swap(ans.Sigma_inv, Sigma_inv);

  if(!is_forward)
    mu = data.err_state->map(mu).sv;
  mu = solve_w_precomputed_chol(ans.Sigma_inv_chol, mu);

  return ans;
}

input_for_normal_apprx taylor_normal_approx(
    pf_base_dens &dens_calc, const PF_data &data, const unsigned int t,
    const arma::mat &Q_inv, const arma::vec &alpha_bar,
    const unsigned int debug_lvl, const bool multithread,
    const bool is_forward){
  /*
   Had similar issues as posted here:
   http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-June/005968.html

   Thus, I made this overload
   */

  arma::uvec r_set = get_risk_set(data, t);

  return(
    taylor_normal_approx
    (dens_calc, data, t, Q_inv, alpha_bar, r_set, debug_lvl, multithread,
     is_forward));
}

/*----------------------------------------------------*/

inline void debug_before_solve(
    const PF_data &data, const arma::vec &parent,
    const arma::mat &A, const arma::vec &x, const arma::vec &xtra_term){
  const int debug_lvl = 5;
  if(data.debug > debug_lvl - 1L){
    data.log(debug_lvl) <<"Finding new particle for" << std::endl
                        << parent.t()
                        << "Solving Ax = b with A" << std::endl
                        << A
                        << "and b^\\top " << std::endl
                        << x.t()
                        << "Yielding" << std::endl
                        << arma::solve(A, x).t()
                        << "Then adding c with c^\\top" << std::endl
                        << xtra_term.t();
  }
}

input_for_normal_apprx_w_cloud_mean
  taylor_normal_approx_w_cloud_mean(
    pf_base_dens &dens_calc, const PF_data &data,
    const unsigned int t, const covarmat &Q, const arma::vec &alpha_bar,
    cloud &cl, const bool is_forward){
    const arma::mat *Q_inv;
    const arma::vec *mu_term;
    if(!is_forward){
      arma::mat tmp = data.state_trans_err->map(Q.inv(), both, trans).sv;
      Q_inv = new arma::mat(tmp + data.uncond_covar_inv(t));
      mu_term = &data.uncond_mean_term(t);

    } else {
      Q_inv = &Q.inv();
      // avoid wmaybe-uninitialized
      mu_term = nullptr;

    }

    input_for_normal_apprx_w_cloud_mean ans =
      taylor_normal_approx
      (dens_calc, data, t, *Q_inv, alpha_bar, 2, true, is_forward);

    auto n_elem = cl.size();
    ans.mu_js = std::vector<arma::vec>(n_elem);
    ans.xi_js = std::vector<arma::vec>(n_elem);
#ifdef _OPENMP
#pragma omp  parallel for schedule(static)
#endif
  for(unsigned int i = 0; i < n_elem; ++i){
    arma::vec &mu_j = ans.mu_js[i];
    arma::vec &xi_j = ans.xi_js[i];
    particle &pr = cl[i];

    if(is_forward){
      mu_j = data.state_trans_err->map(pr.get_state()).sv;
      mu_j = solve_w_precomputed_chol(Q.chol(), mu_j);

    } else {
      mu_j = data.err_state_inv->map(pr.get_state()).sv;
      mu_j = solve_w_precomputed_chol(Q.chol(), mu_j);
      mu_j = data.state_trans_err->map(mu_j, trans).sv + *mu_term;

    }

    debug_before_solve(data, pr.get_state(), ans.Sigma_inv, mu_j, ans.mu);
    mu_j = solve_w_precomputed_chol(ans.Sigma_inv_chol, mu_j) + ans.mu;

    if(is_forward){
      xi_j = std::move(mu_j);
      mu_j = data.state_trans->map(pr.get_state()).sv;
      mu_j.elem(data.err_state->non_zero_row_idx()) =
        xi_j.elem(data.err_state->non_zero_col_idx());

    } else
      xi_j = data.err_state_inv->map(mu_j).sv;

  }

  if(!is_forward){
    delete Q_inv;

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
  const arma::mat *Q_inv;
  arma::mat Q_art_chol;
  const arma::vec *mu_term;
  if(!is_forward){
    // Add the covariance matrix of the artificial prior
    arma::mat tmp = data.state_trans_err->map(Q.inv(), both, trans).sv;
    Q_inv = new arma::mat(tmp + data.uncond_covar_inv(t));
    mu_term = &data.uncond_mean_term(t);

  } else {
    Q_inv = &Q.inv();
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
   arma::vec alpha_bar = this_state;
   if(is_forward)
     alpha_bar = data.state_trans->map(alpha_bar).sv;
   else
     alpha_bar = data.bw_mean(t, alpha_bar);

   auto inter = taylor_normal_approx
     (dens_calc, data, t, *Q_inv, alpha_bar, r_set, 5, false, is_forward);

   arma::vec &mu = ans[i].mu;
   arma::vec &xi = ans[i].xi;
   if(is_forward){
     mu = data.state_trans_err->map(this_state).sv;
     mu = solve_w_precomputed_chol(Q.chol(), mu);

   } else {
     mu = data.err_state_inv   ->map(this_state).sv;
     mu = solve_w_precomputed_chol(Q.chol(), mu);
     mu = data.state_trans_err->map(mu, trans).sv + *mu_term;

   }

   debug_before_solve(data, this_state, inter.Sigma_inv, mu, inter.mu);
   mu = solve_w_precomputed_chol(inter.Sigma_inv_chol, mu) + inter.mu;

   if(is_forward){
     xi = std::move(mu);
     mu = data.state_trans->map(this_state).sv;
     mu.elem(data.err_state->non_zero_row_idx()) =
       xi.elem(data.err_state->non_zero_col_idx());

   } else
     xi = data.err_state_inv->map(mu).sv;

   std::swap(ans[i].sigma_chol_inv, inter.sigma_chol_inv);
   std::swap(ans[i].sigma, inter.Sigma);
  }

  if(!is_forward){
   delete Q_inv;

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
