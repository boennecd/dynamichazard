// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <armadillo>
#include <Rcpp.h>
#include <thread>
#include <future>

#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#define ARMA_HAVE_STD_ISFINITE
#define ARMA_HAVE_STD_ISINF
#define ARMA_HAVE_STD_ISNAN
#define ARMA_HAVE_STD_SNPRINTF


// Rcpp has its own stream object which cooperates more nicely
// with R's i/o -- and as of Armadillo 2.4.3, we can use this
// stream object as well
#if !defined(ARMA_DEFAULT_OSTREAM)
#define ARMA_DEFAULT_OSTREAM Rcpp::Rcout
#endif

// R now defines NDEBUG which suppresses a number of useful
// Armadillo tests Users can still defined it later, and/or
// define ARMA_NO_DEBUG
/*#if defined(NDEBUG)
#undef NDEBUG
#endif*/

// Maybe look at this one day https://github.com/Headtalk/armadillo-ios/blob/master/armadillo-4.200.0/include/armadillo_bits/fn_trunc_exp.hpp
const double lower_trunc_exp_log_thres = sqrt(log(std::numeric_limits<double>::max())) - 1.1;
const double lower_trunc_exp_exp_thres = exp(lower_trunc_exp_log_thres);
inline void in_place_lower_trunc_exp(arma::colvec & result)
{
  for(auto i = result.begin(); i != result.end(); i++){
    if(*i >= lower_trunc_exp_log_thres )
    {
      *i =  lower_trunc_exp_exp_thres;
    }
    else
    {
      *i = std::exp(*i);
    }
  }
}

inline double lower_trunc_exp(double x){
  if(x >= lower_trunc_exp_log_thres )
  {
    return(lower_trunc_exp_exp_thres);
  }
  else
  {
    return(std::exp(x));
  }
}

// from http://gallery.rcpp.org/articles/vector-minimum/
template<typename T>
inline int vecmin(const T x){
  // Rcpp supports STL-style iterators
  auto it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}

// Define convergence criteria
double relative_norm_change(const arma::vec &prev_est, const arma::vec &new_est){
  return arma::norm(prev_est - new_est, 2) / (arma::norm(prev_est, 2) + 1.0e-10);
}
double (*conv_criteria)(const arma::vec&, const arma::vec&) = relative_norm_change;

// locks for parallel implementation
std::mutex m_u;
std::mutex m_U;

//' @export
// [[Rcpp::export]]
Rcpp::List ddhazard_fit_cpp_prelim(const Rcpp::NumericMatrix &X, const arma::vec &tstart,
                                   const arma::vec &tstop, const arma::ivec &events, // armadillo have no boolean vector
                                   const arma::colvec &a_0,
                                   arma::mat Q_0, // by value copy. This  is key cuz we will change it if est_Q_0 = T
                                   arma::mat Q, // similary this is a copy
                                   const Rcpp::List &risk_obj,
                                   const arma::mat &F_,
                                   const int n_max = 100, const double eps = 0.001,
                                   const bool verbose = false, const bool save_all_output = false,
                                   const int order_ = 1, const bool est_Q_0 = true){
  // Initalize constants
  const int d = Rcpp::as<int>(risk_obj["d"]);
  const double Q_warn_eps = sqrt(std::numeric_limits<double>::epsilon());
  const Rcpp::List &risk_sets = Rcpp::as<Rcpp::List>(risk_obj["risk_sets"]);
  const int n_parems = a_0.size() / order_;

  const arma::mat T_F_ = F_.t();
  // Note a copy of data is made below and it is not const for later initalization of pointer to memory (this is not possible with a const pointer)
  arma::mat _X = arma::mat(X.begin(), X.nrow(), X.ncol()).t(); // Armadillo use column major ordering https://en.wikipedia.org/wiki/Row-major_order

  const std::vector<double> I_len = Rcpp::as<std::vector<double> >(risk_obj["I_len"]);

  // Declare and maybe intialize non constants
  double event_time, delta_t, test_max_diff;

  Rcpp::NumericVector conv_values;

  arma::colvec a_prev(a_0.begin(), a_0.size());
  Rcpp::List all_output; // only need if save_all_output = true

  arma::mat a_t_t_s(n_parems * order_, d + 1);
  arma::mat a_t_less_s(n_parems * order_, d);

  arma::cube V_t_t_s(n_parems * order_, n_parems * order_, d + 1);
  arma::cube V_t_less_s(n_parems * order_, n_parems * order_, d);
  arma::cube B_s(n_parems * order_, n_parems * order_, d);

  a_t_t_s.col(0) = a_0;

  arma::colvec u(n_parems * order_);
  arma::mat U(n_parems * order_, n_parems * order_);

  arma::uvec r_set;
  arma::mat V_t_less_s_inv;

  // Needed for lag one covariance
  const int n_in_last_set = (Rcpp::as<arma::uvec>(risk_sets[d - 1])).size();
  arma::mat z_dot(n_parems * order_, n_in_last_set);
  arma::vec H_diag_inv(n_in_last_set);
  arma::mat K_d;
  arma::cube lag_one_cor(n_parems * order_, n_parems * order_, d);

  unsigned int it = 0;

  // M-stp pointers for convenience
  // see http://stackoverflow.com/questions/35983814/access-column-of-matrix-without-creating-copy-of-data-using-armadillo-library
  // the use of unsafe_col is key
  arma::mat *B, *V_less, *V;
  arma::vec a_less, a;

  // Hepler structure to reference data
  struct problem_data{
    const int n_parems;
    arma::mat &_X;
    arma::mat &U;
    arma::vec &u;
    arma::mat &z_dot;
    arma::vec &H_diag_inv;

    const arma::ivec &events;
    const arma::vec &tstop;

    problem_data(arma::mat &X_ref, arma::mat &U_ref,
                 arma::vec &u_ref,
                 arma::mat &z_dot_ref, arma::vec &H_diag_inv_ref,
                 const arma::ivec &events_ref,
                 const arma::vec &tstop_ref,
                 const int n_par):
      _X(X_ref),
      U(U_ref), u(u_ref), n_parems(n_par),
      z_dot(z_dot_ref), H_diag_inv(H_diag_inv_ref),
      events(events_ref), tstop(tstop_ref)
    {}

    const double event_eps = 100 * std::numeric_limits<double>::epsilon(); // something small
  };

  problem_data p_data(_X, U, u, z_dot, H_diag_inv, events, tstop, n_parems); // TODO: use this object through code

  // Class to handle parallel computation of E-step filter step
  class filter_step_helper{
    using uvec_iter = arma::uvec::const_iterator;

    // handy class to ensure that all ressources are cleaned up
    class join_threads
    {
      std::vector<std::thread>& threads;
    public:
      explicit join_threads(std::vector<std::thread>& threads_):
        threads(threads_)
      {}
      ~join_threads()
      {
        for(unsigned long i=0;i<threads.size();++i)
        {
          if(threads[i].joinable())
            threads[i].join();
        }
      }
    };

    // worker class for parallel computation
    class filter_worker{
      bool is_first_call;
      problem_data &dat;

      // local variables to compute temporary result
      arma::colvec u_;
      arma::mat U_;

    public:
      filter_worker(problem_data &p_data):
      dat(p_data), is_first_call(true)
      {}

      bool operator()(uvec_iter first, uvec_iter last,
                    const arma::vec &i_a_t, bool compute_z_and_H,
                    double event_time, int i_start){
        // potentially intialize variables and set entries to zeroes in any case
        if(is_first_call){
          u_ = arma::vec(dat.n_parems);
          U_ = arma::mat(dat.n_parems, dat.n_parems);
          is_first_call = false;
        }
        u_.zeros();
        U_.zeros();

        // compute local results
        int i = i_start;
        for(auto it = first; it != last; it++){
          const arma::vec x_(dat._X.colptr(*it), dat.n_parems, false);
          const double i_eta = lower_trunc_exp(arma::dot(i_a_t, x_));

          if(dat.events(*it) && std::abs(dat.tstop(*it) - event_time) < dat.event_eps){
            u_ += x_ * (1.0 - i_eta / (i_eta + 1.0));
          }
          else {
            u_ -= x_ * (i_eta / (i_eta + 1.0));
          }
          U_ += x_ *  (x_.t() * (i_eta / pow(i_eta + 1.0, 2.0))); // I guess this is the fastest http://stackoverflow.com/questions/26766831/armadillo-inplace-plus-significantly-slower-than-normal-plus-operation

          if(compute_z_and_H){
            dat.H_diag_inv(i) = pow(1.0 + i_eta, 2.0) / i_eta;
            dat.z_dot.rows(0, dat.n_parems - 1).col(i) = x_ *  (i_eta / pow(1.0 + i_eta, 2.0));
            ++i;
          }
        }

        // Update shared variable
        {
          std::lock_guard<std::mutex> lk(m_U);
          dat.U.submat(0, 0, dat.n_parems - 1, dat.n_parems - 1) +=  U_;
        }

        {
          std::lock_guard<std::mutex> lk(m_u);
          dat.u.head(dat.n_parems) += u_;
        }

        return true;
      }
    };

    unsigned long const hardware_threads;
    std::vector<filter_worker> workers;
    problem_data &dat;

  public:
    filter_step_helper(problem_data &p_data):
    hardware_threads(std::thread::hardware_concurrency()),
    dat(p_data), workers()
    {
      // create workers
      for(int i = 0; i < hardware_threads; i++){
        workers.push_back(filter_worker(p_data));
      }
    }

    void parallel_filter_step(uvec_iter first, uvec_iter last,
                              const arma::vec &i_a_t,
                              const bool compute_H_and_z,
                              double event_time){
      // Set referenced objects entries to zero
      dat.U.zeros();
      dat.u.zeros();
      if(compute_H_and_z){
        dat.z_dot.zeros();
        dat.H_diag_inv.zeros();
      }

      // Compute the number of threads to create
      unsigned long const length = std::distance(first, last);

      unsigned long const min_per_thread = 25; // TODO: how to set?
      unsigned long const max_threads =
        (length + min_per_thread - 1) / min_per_thread;

      unsigned long const num_threads =
        std::min(hardware_threads != 0 ? hardware_threads:2, max_threads); // at least use two threads if hardware thread is not set

      unsigned long const block_size = length/num_threads;
      std::vector<std::future<bool> > futures(num_threads - 1);
      std::vector<std::thread> threads(num_threads - 1);
      join_threads joiner(threads);

      // start workers
      // declare outsite for loop to ref after loop
      uvec_iter block_start = first;
      auto it = workers.begin();
      int i_start = 0;
      for(unsigned long i = 0; i < num_threads - 1; ++i, ++it)
      {
        uvec_iter block_end = block_start;
        std::advance(block_end,block_size);
        std::packaged_task<bool(uvec_iter, uvec_iter, const arma::vec&, bool,
                                double, int)> task(*it);
        futures[i] = task.get_future();
        threads[i] = std::thread(std::move(task), block_start, block_end, i_a_t, compute_H_and_z,
                                 event_time, i_start);

        i_start += block_size;
        block_start = block_end;
      }
      (*it)(block_start, last, i_a_t, compute_H_and_z, event_time, i_start); // compute last enteries on this thread

      for(unsigned long i = 0; i < num_threads - 1; ++i)
      {
        futures[i].get(); // will throw if any of the threads did
      }
    }
  };

  //EM algorithm
  filter_step_helper filter_helper(p_data);
  do
  {
    V_t_t_s.slice(0) = Q_0; // Q_0 may have been updated or not

    // E-step
    event_time = vecmin(tstart);

    for (int t = 1; t < d + 1; t++){
      delta_t = I_len[t - 1];
      event_time += delta_t;

      // E-step: Filter step
      a_t_less_s.col(t - 1) = F_ *  a_t_t_s.unsafe_col(t - 1);
      V_t_less_s.slice(t - 1) = F_ * V_t_t_s.slice(t - 1) * T_F_ + delta_t * Q;

      // E-step: scoring step: information matrix and scoring vector
      r_set = Rcpp::as<arma::uvec>(risk_sets[t - 1]) - 1;
      const arma::vec i_a_t(a_t_less_s.colptr(t - 1), n_parems, false);

      filter_helper.parallel_filter_step(r_set.begin(), r_set.end(), i_a_t, t == d, event_time);

      // E-step: scoring step: update values
      V_t_less_s_inv = inv_sympd(V_t_less_s.slice(t - 1));
      V_t_t_s.slice(t) = inv_sympd(V_t_less_s_inv + U);
      a_t_t_s.col(t) = a_t_less_s.unsafe_col(t - 1) + V_t_t_s.slice(t) * u;
      B_s.slice(t - 1) = V_t_t_s.slice(t - 1) * T_F_ * V_t_less_s_inv;

      if(t == d){
        K_d = V_t_less_s.slice(t - 1) * inv(arma::eye<arma::mat>(size(U)) + U * V_t_less_s.slice(t - 1)) * z_dot * diagmat(H_diag_inv);
        // Parenthesis is key here to avoid making a n x n matrix for large n
        K_d = (F_ * V_t_less_s.slice(t - 1) * z_dot * diagmat(H_diag_inv) * z_dot.t()) * K_d;
        K_d = F_ * V_t_less_s.slice(t - 1) * z_dot * diagmat(H_diag_inv) -  K_d;

        lag_one_cor.slice(t - 1) = (arma::eye<arma::mat>(size(U)) - K_d * z_dot.t()) * F_ * V_t_t_s.slice(t - 1);
      }
    }

    // E-step: smoothing
    for (int t = d - 1; t > -1; t--){
      // we need to compute the correlation matrix first
      if(t > 0){
        lag_one_cor.slice(t - 1) = V_t_t_s.slice(t) * B_s.slice(t - 1).t() +  B_s.slice(t) * (
          lag_one_cor.slice(t) - F_ * V_t_t_s.slice(t)) * B_s.slice(t - 1).t();
      }

      a_t_t_s.col(t) = a_t_t_s.unsafe_col(t) + B_s.slice(t) *
        (a_t_t_s.unsafe_col(t + 1) - a_t_less_s.unsafe_col(t));
      V_t_t_s.slice(t) = V_t_t_s.slice(t) + B_s.slice(t) *
        (V_t_t_s.slice(t + 1) - V_t_less_s.slice(t)) * B_s.slice(t).t();

      /*if(t < d + 1 && i_am < 2){ // TODO: Delete
      Rcpp::Rcout << std::setprecision(17) << "It = " << it << " t = " << t << std::endl;
      V_t_t_s.slice(t).raw_print(Rcpp::Rcout);
      a_t_t_s.col(t).raw_print(Rcpp::Rcout);
    }*/
  }

    // M-step
    if(est_Q_0){
      Q_0 = V_t_t_s.slice(0);
    }
    Q.zeros();
    for (int t = 1; t < d + 1; t++){
      delta_t = I_len[t - 1];

      B = &B_s.slice(t - 1);
      V_less = &V_t_t_s.slice(t - 1);
      V = &V_t_t_s.slice(t);
      a_less = a_t_t_s.unsafe_col(t - 1);
      a = a_t_t_s.unsafe_col(t);

      Q += ((a - F_ * a_less) * (a - F_ * a_less).t() + *V
              - F_ * *B * *V
              - (F_ * *B * *V).t()
              + F_ * *V_less * T_F_) / delta_t;
    }
    Q /= d;

    if((test_max_diff = static_cast<arma::mat>(Q - Q.t()).max()) > Q_warn_eps){
      std::ostringstream warning;
      warning << "Q - Q.t() maximal element difference was " << test_max_diff <<
        " in iteration " << it + 1;
      Rcpp::warning(warning.str());
    }

    if((test_max_diff = static_cast<arma::mat>(Q_0 - Q_0.t()).max()) > Q_warn_eps){
      std::ostringstream warning;
      warning << "Q_0 - Q_0.t() maximal element difference was " << test_max_diff <<
        " in iteration " << it + 1;
      Rcpp::warning(warning.str());
    }

    //TODO: better way to ensure that Q and Q_0 are symmetric?
    Q = (Q + Q.t()) / 2.0;
    Q_0 = (Q_0 + Q_0.t()) / 2.0;

    if(order_ > 1){ // CHANGED # TODO: I figure I should set the primaery element to zero, right?
      arma::mat tmp_Q = Q.submat(0, 0, n_parems - 1, n_parems - 1);
      Q.zeros();
      Q.submat(0, 0, n_parems - 1, n_parems - 1) = tmp_Q;
    }

    /*Rcpp::Rcout << std::setprecision(16) << "It " << it + 1 << std::endl; //TODO: delete
    Q.raw_print(Rcpp::Rcout);
    Q_0.raw_print(Rcpp::Rcout);*/

    conv_values.push_back(conv_criteria(a_prev, a_t_t_s.unsafe_col(0)));

    //if(save_all_output) // TODO: make similar save all output function?

    if(*(conv_values.end() -1) < eps){
      break;
    }

    //if(verbose) // TODO: Implement verbose stuff

    a_prev = a_t_t_s.col(0);
}while(++it < n_max);

  /* TODO: Implementt something similar
  // set names
  tmp_names = rep(colnames(X), order_)
  colnames(inner_output$a_t_d_s) = tmp_names
  dimnames(inner_output$V_t_d_s) = list(tmp_names, tmp_names, NULL)
  dimnames(Q) = dimnames(Q_0) = list(tmp_names, tmp_names)
  */

  /* if(save_all_output){
  all_output
  }else */

  if(it == n_max)
    throw std::runtime_error("EM algorithm did not converge within the n_max number of iterations");

  return(Rcpp::List::create(Rcpp::Named("V_t_d_s") = Rcpp::wrap(V_t_t_s),
                            Rcpp::Named("a_t_d_s") = Rcpp::wrap(a_t_t_s.t()),
                            Rcpp::Named("B_s") = Rcpp::wrap(B_s),
                            Rcpp::Named("lag_one_cor") = Rcpp::wrap(lag_one_cor),

                            Rcpp::Named("n_iter") = it + 1,
                            Rcpp::Named("conv_values") = conv_values,
                            Rcpp::Named("Q") = Rcpp::wrap(Q),
                            Rcpp::Named("Q_0") = Rcpp::wrap(Q_0)));
  }


/*** R
getwd()
  library(survival); library(benssurvutils); library(dynamichazard); source("../R/test_utils.R")

    set.seed(2972)
    sims <- test_sim_func_logit(n_series = 10^4, n_vars = 3, t_0 = 0, t_max = 10,
                                x_range = .1, x_mean = -.4, re_draw = T)
    sims$res <- as.data.frame(sims$res)

    design_mat <- benssurvutils::get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res)
    rist_sets <- benssurvutils::get_risk_sets(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

    log_file = file("debug.log")
    sink(log_file)

    tryCatch({
      res_new <- ddhazard_fit_cpp_prelim(
          X = design_mat$X,
          tstart = design_mat$Y[, 1],  tstop = design_mat$Y[, 2], events = design_mat$Y[, 3],
                                                                                       a_0 = rep(0, ncol(design_mat$X)),
                                                                                       Q_0 = diag(10, ncol(design_mat$X)), # something large
                                                                                         Q = diag(1, ncol(design_mat$X)), # something large
                                                                                           F_ = diag(1, ncol(design_mat$X)), # first order random walk
                                                                                             risk_obj = rist_sets,
                                                                                               eps = 10^-4, n_max = 10^4,
                                                                                               order_ = 1,
                                                                                               est_Q_0 = F_)
    }, finally = function(...){
      close(log_file)
      sink()
    })
    */
