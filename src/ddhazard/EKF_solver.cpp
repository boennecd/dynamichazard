#include "../ddhazard.h"
#include "../arma_BLAS_LAPACK.h"
#include "../utils.h"
#include "../family.h"
#include "../thread_pool.h"


class ddhazard_data_EKF {
public:
  ddhazard_data &org;

  // locks for parallel implementation when we need to perform reduction from
  // local score vectors and information matricies
  std::mutex m_u;
  std::mutex m_U;

  const int n_in_last_set;
  const bool is_mult_NR;
  const double NR_eps;
  const unsigned int NR_it_max;

  // Other EKF parameters
  const int EKF_batch_size;

  // Vector for score and information matrix
  arma::colvec u;
  arma::mat U;

  // Needed for lag one covariance
  arma::mat z_dot;
  arma::vec H_diag_inv;
  arma::mat K_d;

  ddhazard_data_EKF(
    ddhazard_data &org, Rcpp::Nullable<Rcpp::NumericVector> NR_eps,
    const unsigned int NR_it_max, const int EKF_batch_size) :
    org(org),
    n_in_last_set(Rcpp::as<arma::uvec>(org.risk_sets[org.d - 1]).size()),
    is_mult_NR(NR_eps.isNotNull()),
    NR_eps(is_mult_NR ? Rcpp::as< Rcpp::NumericVector >(NR_eps)[0] : 0.0),
    NR_it_max(NR_it_max),
    EKF_batch_size(EKF_batch_size){
    if(org.any_dynamic || org.any_fixed_in_E_step){
      u = arma::colvec(org.covar_dim);
      U = arma::mat(org.covar_dim, org.covar_dim);

      z_dot = arma::mat(
        org.state_dim, // TODO: can be changed to dimension of covariates
        n_in_last_set, arma::fill::zeros);
      H_diag_inv = arma::vec(n_in_last_set);
    }
  }
};



class EKF_filter_worker{
  void do_comps(const arma::uvec::const_iterator it, int i,
                const arma::vec &dynamic_coefs, const bool compute_z_and_H,
                const int bin_number,
                const double bin_tstart, const double bin_tstop);

  // Variables for computations
  ddhazard_data_EKF &dat;
  ddhazard_data &org;
  arma::uvec::const_iterator first;
  const arma::uvec::const_iterator last;
  const arma::vec &dynamic_coefs;
  const bool compute_z_and_H;
  const int i_start;
  const int bin_number;
  const double bin_tstart;
  const double bin_tstop;
  family_base &fam;

  // local variables to compute temporary result
  arma::colvec u_;
  arma::mat U_;

public:
  EKF_filter_worker(
    ddhazard_data_EKF &p_data,
    arma::uvec::const_iterator first_, const arma::uvec::const_iterator last_,
    const arma::vec &dynamic_coefs, const bool compute_z_and_H_,
    const int i_start_, const int bin_number_,
    const double bin_tstart_, const double bin_tstop_,
    family_base &fam);

  void operator()();
};

// worker class for parallel computation
EKF_filter_worker::EKF_filter_worker(
  ddhazard_data_EKF &p_data,
  arma::uvec::const_iterator first_, const arma::uvec::const_iterator last_,
  const arma::vec &dynamic_coefs, const bool compute_z_and_H_,
  const int i_start_, const int bin_number_,
  const double bin_tstart_, const double bin_tstop_, family_base &fam)

  :

  dat(p_data), org(p_data.org), first(first_), last(last_),
  dynamic_coefs(dynamic_coefs), compute_z_and_H(compute_z_and_H_),
  i_start(i_start_), bin_number(bin_number_),
  bin_tstart(bin_tstart_), bin_tstop(bin_tstop_), fam(fam),
  u_(org.covar_dim, arma::fill::zeros),
  U_(org.covar_dim, org.covar_dim, arma::fill::zeros)
{}

inline void EKF_filter_worker::operator()(){
  // potentially intialize variables and set entries to zeroes in any case

  // compute local results
  int i = i_start;
  const bool uses_at_risk_length = fam.uses_at_risk_length();
  for(arma::uvec::const_iterator it = first; it != last; it++, i++){
    const arma::vec x_(org.X.colptr(*it), org.covar_dim, false);
    const double w = org.weights(*it);
    const double offset = org.fixed_effects(*it);
    const double eta = arma::dot(dynamic_coefs, x_) + offset;
    const bool do_die = org.is_event_in_bin(*it) == bin_number;
    const double at_risk_length =
      uses_at_risk_length ?
        get_at_risk_length(
          org.tstop(*it), bin_tstop, org.tstart(*it), bin_tstart) :
        0.;

    auto trunc_res = fam.truncate_eta(do_die, eta, exp(eta), at_risk_length);
    const double mu = fam.linkinv(trunc_res, at_risk_length);
    const double mu_eta = fam.mu_eta(trunc_res, at_risk_length);
    const double var = fam.var(trunc_res, at_risk_length);

    /* Update local score and Hessian */
    u_ += x_ * (w * mu_eta * (do_die - mu) / (var + org.denom_term));
    sym_mat_rank_one_update(
      w * mu_eta * mu_eta / (var + org.denom_term), x_, U_);

    if(compute_z_and_H){
      dat.H_diag_inv(i) = 1 / var;
      arma::vec tmp(x_ * mu_eta);
      dat.z_dot.col(i) = org.state_lp_inv->map(tmp).sv;
    }
  }

  // Update shared variable
  {
    std::lock_guard<std::mutex> lk(dat.m_U);
    dat.U += U_;
  }

  {
    std::lock_guard<std::mutex> lk(dat.m_u);
    dat.u += u_;
  }
}


EKF_solver::EKF_solver(
  ddhazard_data &p, const std::string model,
  Rcpp::Nullable<Rcpp::NumericVector> NR_eps,
  const unsigned int NR_it_max, const int EKF_batch_size, family_base &fam) :
  org(p), p_dat(new ddhazard_data_EKF(
      p, NR_eps, NR_it_max, EKF_batch_size)), pool(new thread_pool(org.n_threads)),
      model(model), fam(fam)
  {}

void EKF_solver::solve(){
  double bin_tstop = org.min_start;

  for (int t = 1; t < org.d + 1; t++){

    double bin_tstart = bin_tstop;
    double delta_t = org.I_len[t - 1];
    bin_tstop += delta_t;
    // E-step: Prediction step
    org.a_t_less_s.col(t - 1) =
      org.state_trans->map(org.a_t_t_s.col(t - 1)).sv;

    org.V_t_less_s.slice(t - 1) =
      org.state_trans->map(org.V_t_t_s.slice(t - 1)).sv +
      delta_t * org.err_state->map(org.Q).sv;

    if(org.debug){
      std::stringstream str;
      str << t << "|" << t - 1;

      my_print(org, org.a_t_less_s.col(t - 1), "a_(" + str.str() + ")");
      my_print(org, org.V_t_less_s.slice(t - 1), "V_(" + str.str() + ")");
      my_debug_logger(org)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(org.V_t_less_s.slice(t - 1));
    }

    // E-step: scoring step: information matrix and scoring vector
    arma::uvec r_set = get_risk_set(org, t);
    arma::vec i_a_t = org.a_t_less_s.col(t - 1);
    arma::mat V_t_less_s_inv;
    inv_sympd(V_t_less_s_inv, org.V_t_less_s.slice(t - 1), org.use_pinv,
              "ddhazard_fit_cpp estimation error: Failed to invert V_(t|t-1)");
    unsigned int n_NR_it = 0;

    arma::vec update_term = V_t_less_s_inv * org.a_t_less_s.col(t - 1);

    while(true){
      ++n_NR_it;

      arma::vec dynamic_coef(org.state_lp->map(i_a_t).sv);
      parallel_filter_step(
        r_set.begin(), r_set.end(), dynamic_coef,
        t == org.d, t - 1, bin_tstart, bin_tstop);

      if(org.debug){
        my_debug_logger(org) << "Score vector and information matrix at time " << t << " are:";
        my_print(org, p_dat->u, "u");
        my_print(org, p_dat->U, "U");
      }

      if(p_dat->u.has_inf() || p_dat->u.has_nan()){
        Rcpp::stop("ddhazard_fit_cpp estimation error: score vector had inf or nan elements. Try decreasing the learning rate");

      } else if(p_dat->U.has_inf() || p_dat->U.has_nan()){
        Rcpp::stop("ddhazard_fit_cpp estimation error: information matrix had inf or nan elements. Try decreasing the learning rate");

      }

      // E-step: scoring step: update values
      auto U = org.state_lp_inv->map(p_dat->U);
      inv_sympd(
        org.V_t_t_s.slice(t) , V_t_less_s_inv + U.sv,
        org.use_pinv,
        "ddhazard_fit_cpp estimation error: Failed to compute inverse for V_(t|t)");

      auto u = org.state_lp_inv->map(p_dat->u);
      org.a_t_t_s.col(t) =
        org.V_t_t_s.slice(t) * (
            U.sv * i_a_t + update_term + (org.LR * u.sv));

      if(org.debug){
        my_print(org, i_a_t, "a^(" + std::to_string(n_NR_it - 1L) + ")");
        my_print(org, org.a_t_t_s.col(t),
                 "a^(" +  std::to_string(n_NR_it) + ")");
      }

      if(!p_dat->is_mult_NR || arma::norm(org.a_t_t_s.col(t) - i_a_t, 2) /
         (arma::norm(i_a_t, 2) + 1e-8) < p_dat->NR_eps)
        break;

      if(n_NR_it > p_dat->NR_it_max)
        Rcpp::stop("Failed to convergece in NR method of filter step within " +
          std::to_string(p_dat->NR_it_max) + " iterations");

      if(org.debug){
        my_debug_logger(org)
          << "Did not converge in filter step in iteration " << n_NR_it <<
        ". Convergence criteria value is  "
          << arma::norm(org.a_t_t_s.col(t) - i_a_t, 2) /
          (arma::norm(i_a_t, 2) + 1e-8);
      }

      i_a_t = org.a_t_t_s.col(t);
    }

    org.B_s.slice(t - 1) =
      org.state_trans->map(org.V_t_t_s.slice(t - 1), right).sv * V_t_less_s_inv;

    if(org.debug){
      std::stringstream str;
      str << t << "|" << t;

      my_print(org, org.a_t_t_s.col(t), "a_(" + str.str() + ")");
      my_print(org, org.V_t_t_s.slice(t), "V_(" + str.str() + ")");
      my_debug_logger(org)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(org.V_t_t_s.slice(t));
    }

    if(t == org.d){
      auto U = org.state_lp_inv->map(p_dat->U);
      arma::mat tmp_inv_mat;
      inv(tmp_inv_mat,
          arma::eye<arma::mat>(org.state_dim, org.state_dim) +
            U.sv * org.V_t_less_s.slice(t - 1),
          org.use_pinv, "ddhazard_fit_cpp estimation error: Failed to invert intermediate for K_d matrix");

      p_dat->K_d = org.V_t_less_s.slice(t - 1) * (
        tmp_inv_mat * p_dat->z_dot * diagmat(p_dat->H_diag_inv));
      // Parenthesis is key here to avoid making a n x n matrix for large n
      p_dat->K_d = (org.state_trans->map(
        org.V_t_less_s.slice(t - 1), left).sv *
        p_dat->z_dot  * diagmat(p_dat->H_diag_inv) *
        p_dat->z_dot.t()) * p_dat->K_d;
      p_dat->K_d =
        org.state_trans->map(org.V_t_less_s.slice(t - 1), left).sv *
        p_dat->z_dot  * diagmat(p_dat->H_diag_inv) -  p_dat->K_d;

      org.lag_one_cov.slice(t - 1) =
        (arma::eye<arma::mat>(org.state_dim, org.state_dim) -
        p_dat->K_d * p_dat->z_dot.t()) *
        org.state_trans->map(org.V_t_t_s.slice(t - 1), left).sv;
    }
  }
}

void EKF_solver::parallel_filter_step(
    arma::uvec::const_iterator first, arma::uvec::const_iterator last,
    const arma::vec &dynamic_coefs,
    const bool compute_H_and_z,
    const int bin_number,
    const double bin_tstart, const double bin_tstop){

  // Set entries to zero
  p_dat->U.zeros();
  p_dat->u.zeros();
  if(compute_H_and_z){
    p_dat->z_dot.zeros();
    p_dat->H_diag_inv.zeros();
  }

  // Compute the number of blocks to create
  unsigned long const length = std::distance(first, last);

  unsigned long const block_size =
    org.n_threads <= 1 ?
    length :
    std::max(p_dat->EKF_batch_size, (int)std::ceil((double)length / org.n_threads));
  unsigned long const num_blocks= (int)std::ceil((double)length / block_size);
  std::vector<std::future<void> > futures(num_blocks - 1);

  std::vector<EKF_filter_worker> workers;
  workers.reserve(num_blocks - 1);

  // declare outsite of loop to ref after loop
  arma::uvec::const_iterator block_start = first;
  int i_start = 0;

  for(unsigned long i = 0; i < num_blocks - 1; ++i){
    arma::uvec::const_iterator block_end = block_start;
    std::advance(block_end, block_size);

    workers.emplace_back(
      *p_dat.get(), block_start, block_end, dynamic_coefs, compute_H_and_z,
      i_start, bin_number, bin_tstart, bin_tstop, fam);

    futures[i] = pool->submit(workers.back());
    i_start += block_size;
    block_start = block_end;
  }

  EKF_filter_worker( // compute last enteries on this thread
    *p_dat.get(), block_start, last, dynamic_coefs, compute_H_and_z,
    i_start, bin_number, bin_tstart, bin_tstop, fam)();

  for(unsigned long i = 0; i < num_blocks - 1; ++i)
  {
    futures[i].get();   // will throw if any of the threads did
  }

  // reflecting the upper triangle to the lower triangle as we have used the
  // dsyr BLAS function
  p_dat->U = symmatu(p_dat->U);
}
