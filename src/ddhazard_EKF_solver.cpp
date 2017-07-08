#include "ddhazard.h"
#include "exp_model_funcs.h"
#include "thread_pool.h"
#include "arma_utils.h"

// worker class for parallel computation
// This class is abstact as the method do_computation will differ between
// the models
EKF_filter_worker::EKF_filter_worker(
  problem_data_EKF &p_data,
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

void EKF_filter_worker::operator()(){
  // potentially intialize variables and set entries to zeroes in any case

  // compute local results
  int i = i_start;
  for(arma::uvec::const_iterator it = first; it != last; it++, i++){
    do_comps(it, i, i_a_t, compute_z_and_H, bin_number,
             bin_tstart, bin_tstop);
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

// worker for the logit model
class EKF_filter_worker_logit : public EKF_filter_worker {
private:
  void do_comps(const arma::uvec::const_iterator it, int i,
                const arma::vec &i_a_t, const bool compute_z_and_H,
                const int bin_number,
                const double bin_tstart, const double bin_tstop){
    const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);
    const double w = dat.weights(*it);

    double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
    const double exp_eta = exp(arma::dot(i_a_t, x_) + offset);

    // Can be issue here with overflow in denominator
    double tmp_denom = pow(1.0 + exp_eta, 2.0);
    const double var = std::isinf(tmp_denom) ?
    pow(exp_eta, -1) : (exp_eta / pow(exp_eta + 1.0, 2.0));

    double score_fac = exp_eta/(dat.denom_term + exp_eta +
                                2*dat.denom_term*exp_eta + dat.denom_term*pow(exp_eta,2));
    double var_fac = pow(exp_eta,2)/(pow(1 + exp_eta,4)*(dat.denom_term + exp_eta/pow(1 + exp_eta, 2)));

    u_ += x_ * (w * score_fac *
      ((dat.is_event_in_bin(*it) == bin_number) - exp_eta / (1.0 + exp_eta)));
    sym_mat_rank_one_update(w * var_fac, x_, U_);

    if(compute_z_and_H){
      dat.H_diag_inv(i) = pow(var, -1);
      dat.z_dot(*dat.span_current_cov, i) = x_ *  var;
    }
  }

public:
  using EKF_filter_worker::EKF_filter_worker;
};


// worker for the continous model with exponential distribution where only the
// binary variable is used
class EKF_filter_worker_exp_bin : public EKF_filter_worker {
private:
  void do_comps(const arma::uvec::const_iterator it, int i,
                const arma::vec &i_a_t, const bool compute_z_and_H,
                const int bin_number,
                const double bin_tstart, const double bin_tstop){
    // Compute intermediates
    const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);
    const double w = dat.weights(*it);

    double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
    const double eta = arma::dot(i_a_t, x_) + offset;

    const double do_die = (dat.is_event_in_bin(*it) == bin_number);
    const double at_risk_length = std::min(dat.tstop(*it), bin_tstop) - std::max(dat.tstart(*it), bin_tstart);

    const double exp_eta = exp(eta);
    const double v = at_risk_length * exp_eta;
    const double exp_v = exp(v);
    const double inv_exp_v = pow(exp_v, -1.0);

    const double expect_chance_die = exp_model_funcs::expect_chance_die(v, inv_exp_v);

    const double score_fac = exp_model_funcs::binary_score_fac(v, inv_exp_v, dat.denom_term);
    const double var_fac = exp_model_funcs::binary_var_fac(v, inv_exp_v, dat.denom_term);

    u_ += x_ * (w * score_fac * (do_die - expect_chance_die));

    sym_mat_rank_one_update(w * var_fac, x_, U_);

    if(compute_z_and_H){
      // Compute terms from waiting time
      dat.H_diag_inv(i) = exp_model_funcs::inv_var_chance_die(v, expect_chance_die);

      dat.z_dot(*dat.span_current_cov, i) =  x_ * (inv_exp_v * v);
    }
  }

public:
  using EKF_filter_worker::EKF_filter_worker;
};

// worker for the continous model with exponential distribution where only the
// right clipped variable is used
class EKF_filter_worker_exp_clip_time : public EKF_filter_worker {
private:
  void do_comps(const arma::uvec::const_iterator it, int i,
                const arma::vec &i_a_t, const bool compute_z_and_H,
                const int bin_number,
                const double bin_tstart, const double bin_tstop){
    // Compute intermediates
    const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);
    const double w = dat.weights(*it);

    double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
    const double eta = arma::dot(i_a_t, x_) + offset;

    const double do_die = (dat.is_event_in_bin(*it) == bin_number);
    const double time_outcome = std::min(dat.tstop(*it), bin_tstop) - std::max(dat.tstart(*it), bin_tstart);
    const double at_risk_length = do_die ? bin_tstop - std::max(dat.tstart(*it), bin_tstart) : time_outcome;

    const double exp_eta = exp(eta);
    const double inv_exp_eta = pow(exp_eta, -1);
    const double v = at_risk_length * exp_eta;
    const double exp_v = exp(v);
    const double inv_exp_v = pow(exp_v, -1.0);

    const double expect_time = exp_model_funcs::expect_time(v, at_risk_length, inv_exp_v, exp_eta);

    const double score_fac = exp_model_funcs::clip_time_score_fac(
      exp_eta, inv_exp_eta, inv_exp_v, at_risk_length, dat.denom_term);
    const double var_fac = exp_model_funcs::clip_time_var_fac(
      exp_eta, inv_exp_eta, inv_exp_v, at_risk_length, dat.denom_term);


    u_ += x_ * (w * score_fac * (time_outcome - expect_time));

    sym_mat_rank_one_update(w * var_fac, x_, U_);

    if(compute_z_and_H){
      // Compute terms from waiting time
      dat.H_diag_inv(i) = exp_model_funcs::inv_var_wait_time(v, exp_eta, inv_exp_v);

      dat.z_dot(*dat.span_current_cov, i) =  x_ * (inv_exp_eta*inv_exp_v*(1 - exp_v + v));
    }
  }

public:
  using EKF_filter_worker::EKF_filter_worker;
};

// worker for the continous model with exponential distribution where only the
// right clipped variable is used where outcomes are negative
class EKF_filter_worker_exp_clip_time_w_jump : public EKF_filter_worker {
private:
  void do_comps(const arma::uvec::const_iterator it, int i,
                const arma::vec &i_a_t, const bool compute_z_and_H,
                const int bin_number,
                const double bin_tstart, const double bin_tstop){
    // Compute intermediates
    const arma::vec x_(dat.X.colptr(*it), dat.n_params_state_vec, false);
    const double w = dat.weights(*it);

    double offset = (dat.any_fixed_in_M_step) ? arma::dot(dat.fixed_parems, dat.fixed_terms.col(*it)) : 0.;
    const double eta = arma::dot(i_a_t, x_) + offset;

    const double do_die = (dat.is_event_in_bin(*it) == bin_number);
    const double time_outcome = std::min(dat.tstop(*it), bin_tstop) - std::max(dat.tstart(*it), bin_tstart);
    const double at_risk_length = do_die ? bin_tstop - std::max(dat.tstart(*it), bin_tstart) : time_outcome;

    const double exp_eta = exp(eta);
    const double inv_exp_eta = pow(exp_eta, -1);
    const double v = at_risk_length * exp_eta;
    const double exp_v = exp(v);
    const double inv_exp_v = pow(exp_v, -1);

    const double expect_time = exp_model_funcs::expect_time_w_jump(exp_eta, inv_exp_eta, inv_exp_v, at_risk_length);

    const double score_fac = exp_model_funcs::clip_time_w_jump_score_fac(
      exp_eta, v, inv_exp_v, at_risk_length, dat.denom_term);
    const double var_fac = exp_model_funcs::clip_time_w_jump_var_fac(
      exp_eta, v, inv_exp_v, at_risk_length, dat.denom_term);


    u_ += x_ * (w * score_fac * (time_outcome - expect_time - at_risk_length * do_die));

    sym_mat_rank_one_update(w * var_fac, x_, U_);

    if(compute_z_and_H){
      // Compute terms from waiting time
      dat.H_diag_inv(i) =
        1 / exp_model_funcs::var_wait_time_w_jump(exp_eta, inv_exp_v, inv_exp_v);

      dat.z_dot(*dat.span_current_cov, i) =  x_ * (
        -inv_exp_eta + at_risk_length*inv_exp_v - pow(at_risk_length,2)*exp_eta*inv_exp_v + inv_exp_eta*inv_exp_v
      );
    }
  }

public:
  using EKF_filter_worker::EKF_filter_worker;
};

















EKF_solver::EKF_solver(problem_data_EKF &p_, const std::string model_):
  p_dat(p_), model(model_)
  {};

void EKF_solver::solve(){
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
    arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;
    arma::vec i_a_t = p_dat.a_t_less_s.col(t - 1);
    arma::mat V_t_less_s_inv;
    inv_sympd(V_t_less_s_inv, p_dat.V_t_less_s.slice(t - 1), p_dat.use_pinv,
              "ddhazard_fit_cpp estimation error: Failed to invert V_(t|t-1)");
    unsigned int n_NR_it = 0;

    arma::vec update_term = V_t_less_s_inv * p_dat.a_t_less_s.col(t - 1);

    while(true){
      ++n_NR_it;

#if defined(USE_OPEN_BLAS)
      openblas_set_num_threads(1);
#endif

      std::unique_ptr<EKF_helper_base> helper;
      if(model == "logit"){
        helper.reset(new EKF_helper<EKF_filter_worker_logit>(p_dat));

      } else if (model == "exp_bin"){
        helper.reset(new EKF_helper<EKF_filter_worker_exp_bin>(p_dat));

      } else if(model == "exp_clip_time"){
        helper.reset(new EKF_helper<EKF_filter_worker_exp_clip_time>(p_dat));

      } else if(model == "exp_clip_time_w_jump"){
        helper.reset(new EKF_helper<EKF_filter_worker_exp_clip_time_w_jump>(p_dat));

      } else
        Rcpp::stop("EKF is not implemented for model '" + model  +"'");

      helper->parallel_filter_step(
        r_set.begin(), r_set.end(), i_a_t(*p_dat.span_current_cov),
        t == p_dat.d, t - 1, bin_tstart, bin_tstop);

      if(p_dat.debug){
        my_debug_logger(p_dat) << "Score vector and diagonal of information matrix at time " << t << " are:";
        my_print(p_dat, p_dat.u, "u");
        my_print(p_dat, p_dat.U, "U");
      }

#ifdef USE_OPEN_BLAS
      openblas_set_num_threads(p_dat.n_threads);
#endif

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




template<class T>
EKF_helper<T>::EKF_helper(problem_data_EKF &p_data_):
  max_threads((p_data_.n_threads > 1) ? p_data_.n_threads - 1 : 1),
  p_data(p_data_)
{
  if(p_data.debug){
    my_debug_logger(p_data) << "EKF solver will use at most " << max_threads << " threads";
  }
}

template<class T>
void EKF_helper<T>::parallel_filter_step(
    arma::uvec::const_iterator first, arma::uvec::const_iterator last,
    const arma::vec &i_a_t,
    const bool compute_H_and_z,
    const int bin_number,
    const double bin_tstart, const double bin_tstop){

  // Set entries to zero
  p_data.U.zeros();
  p_data.u.zeros();
  if(compute_H_and_z){
    p_data.z_dot.zeros();
    p_data.H_diag_inv.zeros();
  }

  // Compute the number of blocks to create
  unsigned long const length = std::distance(first, last);

  unsigned long const block_size = p_data.EKF_batch_size;
  unsigned long const num_blocks=(length+block_size-1)/block_size;
  std::vector<std::future<void> > futures(num_blocks-1);
  thread_pool pool(num_blocks - 1, max_threads);

  std::vector<T> workers;

  // declare outsite of loop to ref after loop
  arma::uvec::const_iterator block_start = first;
  auto it = workers.begin();
  int i_start = 0;

  for(unsigned long i = 0; i < num_blocks - 1; ++i, ++it){
    arma::uvec::const_iterator block_end = block_start;
    std::advance(block_end, block_size);

    workers.emplace_back(
      p_data, block_start, block_end, i_a_t, compute_H_and_z,
      i_start, bin_number, bin_tstart, bin_tstop);

    futures[i] = pool.submit(workers.back());
    i_start += block_size;
    block_start = block_end;
  }

  T( // compute last enteries on this thread
    p_data, block_start, last, i_a_t, compute_H_and_z,
    i_start, bin_number, bin_tstart, bin_tstop)();

  for(unsigned long i = 0; i < num_blocks - 1; ++i)
  {
    futures[i].get();   // will throw if any of the threads did
  }

  // reflecting the upper triangle to the lower triangle as we have used the
  // dsyr BLAS function
  p_data.U = symmatu(p_data.U);
};








