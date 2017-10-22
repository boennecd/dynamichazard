#include "../ddhazard.h"
#include "../exp_model_funcs.h"
#include "../arma_BLAS_LAPACK.h"
#include "../utils.h"
#include "../family.h"
#include "../thread_pool.h"

// worker class for parallel computation
template<typename T>
EKF_filter_worker<T>::EKF_filter_worker(
  ddhazard_data_EKF &p_data,
  arma::uvec::const_iterator first_, const arma::uvec::const_iterator last_,
  const arma::vec &i_a_t_, const bool compute_z_and_H_,
  const int i_start_, const int bin_number_,
  const double bin_tstart_, const double bin_tstop_)

  :

  dat(p_data), first(first_), last(last_),
  i_a_t(i_a_t_), compute_z_and_H(compute_z_and_H_),
  i_start(i_start_), bin_number(bin_number_),
  bin_tstart(bin_tstart_), bin_tstop(bin_tstop_),
  u_(dat.n_params_state_vec, arma::fill::zeros),
  U_(dat.n_params_state_vec, dat.n_params_state_vec, arma::fill::zeros)
{};

template<typename T>
inline void EKF_filter_worker<T>::operator()(){
  // potentially intialize variables and set entries to zeroes in any case

  // compute local results
  int i = i_start;
  for(arma::uvec::const_iterator it = first; it != last; it++, i++){
    const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);
    const double w = dat.weights(*it);
    const double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
    const double eta = arma::dot(i_a_t, x_) + offset;
    const bool do_die = dat.is_event_in_bin(*it) == bin_number;
    const double at_risk_length =
      T::uses_at_risk_length ?
        get_at_risk_length(
          dat.tstop(*it), bin_tstop, dat.tstart(*it), bin_tstart) :
        0.;

    auto trunc_res = T::truncate_eta(do_die, eta, exp(eta), at_risk_length);
    const double mu = T::linkinv(trunc_res, at_risk_length);
    const double mu_eta = T::mu_eta(trunc_res, at_risk_length);
    const double var = T::var(trunc_res, at_risk_length);

    /* Update local score and Hessian */
    u_ += x_ * (w * mu_eta * (do_die - mu) / (var + dat.denom_term));
    sym_mat_rank_one_update(
      w * mu_eta * mu_eta / (var + dat.denom_term), x_, U_);

    if(compute_z_and_H){
      dat.H_diag_inv(i) = 1 / var;
      dat.z_dot(*dat.span_current_cov, i) = x_ * mu_eta;
    }
  }

  // Update shared variable
  {
    std::lock_guard<std::mutex> lk(dat.m_U);
    dat.U(*dat.span_current_cov, *dat.span_current_cov) +=  U_;
  }

  {
    std::lock_guard<std::mutex> lk(dat.m_u);
    dat.u(*dat.span_current_cov) += u_;
  }
};




template<typename T>
EKF_solver<T>::EKF_solver(ddhazard_data_EKF &p_, const std::string model_):
  p_dat(p_), model(model_),
  max_threads((p_.n_threads > 1) ? p_.n_threads - 1 : 1)
  {};

template<typename T>
void EKF_solver<T>::solve(){
  double bin_tstop = p_dat.min_start;

  for (int t = 1; t < p_dat.d + 1; t++){

    double bin_tstart = bin_tstop;
    double delta_t = p_dat.I_len[t - 1];
    bin_tstop += delta_t;

    // E-step: Prediction step
    p_dat.a_t_less_s.col(t - 1) = p_dat.F_ *  p_dat.a_t_t_s.unsafe_col(t - 1);
    p_dat.V_t_less_s.slice(t - 1) = p_dat.F_ * p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ + delta_t * p_dat.Q;

    if(p_dat.debug){
      std::stringstream str;
      str << t << "|" << t - 1;

      my_print(p_dat, p_dat.a_t_less_s.col(t - 1), "a_(" + str.str() + ")");
      my_print(p_dat, p_dat.V_t_less_s.slice(t - 1), "V_(" + str.str() + ")");
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(p_dat.V_t_less_s.slice(t - 1));
    }

    // E-step: scoring step: information matrix and scoring vector
    arma::uvec r_set = get_risk_set(p_dat, t);
    arma::vec i_a_t = p_dat.a_t_less_s.col(t - 1);
    arma::mat V_t_less_s_inv;
    inv_sympd(V_t_less_s_inv, p_dat.V_t_less_s.slice(t - 1), p_dat.use_pinv,
              "ddhazard_fit_cpp estimation error: Failed to invert V_(t|t-1)");
    unsigned int n_NR_it = 0;

    arma::vec update_term = V_t_less_s_inv * p_dat.a_t_less_s.col(t - 1);

    while(true){
      ++n_NR_it;

      parallel_filter_step(
        r_set.begin(), r_set.end(), i_a_t(*p_dat.span_current_cov),
        t == p_dat.d, t - 1, bin_tstart, bin_tstop);

      if(p_dat.debug){
        my_debug_logger(p_dat) << "Score vector and diagonal of information matrix at time " << t << " are:";
        my_print(p_dat, p_dat.u, "u");
        my_print(p_dat, p_dat.U, "U");
      }

      if(p_dat.u.has_inf() || p_dat.u.has_nan()){
        Rcpp::stop("ddhazard_fit_cpp estimation error: score vector had inf or nan elements. Try decreasing the learning rate");

      } else if(p_dat.U.has_inf() || p_dat.U.has_nan()){
        Rcpp::stop("ddhazard_fit_cpp estimation error: information matrix had inf or nan elements. Try decreasing the learning rate");

      }

      // E-step: scoring step: update values
      inv_sympd(p_dat.V_t_t_s.slice(t) , V_t_less_s_inv + p_dat.U, p_dat.use_pinv,
                "ddhazard_fit_cpp estimation error: Failed to compute inverse for V_(t|t)");

      //p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.col(t - 1) + p_dat.LR * p_dat.V_t_t_s.slice(t) * p_dat.u;
      p_dat.a_t_t_s.col(t) = p_dat.V_t_t_s.slice(t) * (
        p_dat.U * i_a_t + update_term + (p_dat.LR * p_dat.u));

      if(p_dat.debug){
        my_print(p_dat,i_a_t, "a^(" + std::to_string(n_NR_it - 1L) + ")");
        my_print(p_dat, p_dat.a_t_t_s.col(t), "a^(" +  std::to_string(n_NR_it) + ")");
      }

      if(!p_dat.is_mult_NR || arma::norm(p_dat.a_t_t_s.col(t) - i_a_t, 2) / (arma::norm(i_a_t, 2) + 1e-8) < p_dat.NR_eps)
        break;

      if(n_NR_it > p_dat.NR_it_max)
        Rcpp::stop("Failed to convergece in NR method of filter step within " + std::to_string(p_dat.NR_it_max) + " iterations");

      if(p_dat.debug){
        my_debug_logger(p_dat)
          << "Did not converge in filter step in iteration " << n_NR_it << ". Convergence criteria value is  "
          << arma::norm(p_dat.a_t_t_s.col(t) - i_a_t, 2) / (arma::norm(i_a_t, 2) + 1e-8);
      }

      i_a_t = p_dat.a_t_t_s.col(t);
    }

    p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * V_t_less_s_inv;

    if(p_dat.debug){
      std::stringstream str;
      str << t << "|" << t;

      my_print(p_dat, p_dat.a_t_t_s.col(t), "a_(" + str.str() + ")");
      my_print(p_dat, p_dat.V_t_t_s.slice(t), "V_(" + str.str() + ")");
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(p_dat.V_t_t_s.slice(t));
    }

    if(t == p_dat.d){
      arma::mat tmp_inv_mat;
      inv(tmp_inv_mat, arma::eye<arma::mat>(size(p_dat.U)) + p_dat.U * p_dat.V_t_less_s.slice(t - 1),
          p_dat.use_pinv, "ddhazard_fit_cpp estimation error: Failed to invert intermediate for K_d matrix");

      p_dat.K_d = p_dat.V_t_less_s.slice(t - 1) * (tmp_inv_mat * p_dat.z_dot * diagmat(p_dat.H_diag_inv));
      // Parenthesis is key here to avoid making a n x n matrix for large n
      p_dat.K_d = (p_dat.F_ * p_dat.V_t_less_s.slice(t - 1) * p_dat.z_dot  * diagmat(p_dat.H_diag_inv) * p_dat.z_dot.t()) * p_dat.K_d;
      p_dat.K_d = p_dat.F_ * p_dat.V_t_less_s.slice(t - 1) * p_dat.z_dot  * diagmat(p_dat.H_diag_inv) -  p_dat.K_d;

      p_dat.lag_one_cov.slice(t - 1) = (arma::eye<arma::mat>(size(p_dat.U)) - p_dat.K_d * p_dat.z_dot.t()) * p_dat.F_ * p_dat.V_t_t_s.slice(t - 1);
    }
  }
}

template<typename T>
void EKF_solver<T>::parallel_filter_step(
    arma::uvec::const_iterator first, arma::uvec::const_iterator last,
    const arma::vec &i_a_t,
    const bool compute_H_and_z,
    const int bin_number,
    const double bin_tstart, const double bin_tstop){
  using worker_T = EKF_filter_worker<T>;

  // Set entries to zero
  p_dat.U.zeros();
  p_dat.u.zeros();
  if(compute_H_and_z){
    p_dat.z_dot.zeros();
    p_dat.H_diag_inv.zeros();
  }

  // Compute the number of blocks to create
  unsigned long const length = std::distance(first, last);

  unsigned long const block_size =
    p_dat.n_threads <= 1 ?
    length :
    std::max(p_dat.EKF_batch_size, (int)std::ceil((double)length / p_dat.n_threads));
  unsigned long const num_blocks= (int)std::ceil((double)length / block_size);
  std::vector<std::future<void> > futures(num_blocks - 1);
  thread_pool pool(num_blocks-1);

  std::vector<worker_T> workers;

  // declare outsite of loop to ref after loop
  arma::uvec::const_iterator block_start = first;
  int i_start = 0;

  for(unsigned long i = 0; i < num_blocks - 1; ++i){
    arma::uvec::const_iterator block_end = block_start;
    std::advance(block_end, block_size);

    workers.emplace_back(
      p_dat, block_start, block_end, i_a_t, compute_H_and_z,
      i_start, bin_number, bin_tstart, bin_tstop);

    futures[i] = pool.submit(workers.back());
    i_start += block_size;
    block_start = block_end;
  }

  worker_T( // compute last enteries on this thread
    p_dat, block_start, last, i_a_t, compute_H_and_z,
    i_start, bin_number, bin_tstart, bin_tstop)();

  for(unsigned long i = 0; i < num_blocks - 1; ++i)
  {
    futures[i].get();   // will throw if any of the threads did
  }

  // reflecting the upper triangle to the lower triangle as we have used the
  // dsyr BLAS function
  p_dat.U = symmatu(p_dat.U);
};

// Define classes
template class EKF_solver<logistic>;
template class EKF_solver<exponential>;






