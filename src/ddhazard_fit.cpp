// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <omp.h>
#include <iostream>
#include <armadillo>
#include <Rcpp.h>

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
  const double event_eps = d * std::numeric_limits<double>::epsilon();
  const double Q_warn_eps = sqrt(std::numeric_limits<double>::epsilon());
  const Rcpp::List &risk_sets = Rcpp::as<Rcpp::List>(risk_obj["risk_sets"]);
  const int n_parems = a_0.size() / order_;

  const arma::mat T_F_ = F_.t();
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
  arma::mat z_dot;
  arma::vec H_diag_inv;
  arma::mat K_d;
  arma::cube lag_one_cor(n_parems * order_, n_parems * order_, d);

  // Parallel settings and variables
  const int n_threads = omp_get_num_procs() - 1;
  omp_set_num_threads(n_threads);
  int i_am, i_points, i_start, n_cols, n_threads_current;

  unsigned int it = 0;
  arma::colvec u_(n_parems);
  arma::mat U_(n_parems, n_parems);
  bool is_run_parallel;

  // M-stp pointers for convenience
  // see http://stackoverflow.com/questions/35983814/access-column-of-matrix-without-creating-copy-of-data-using-armadillo-library
  // the use of unsafe_col is key
  arma::mat *B, *V_less, *V;
  arma::vec a_less, a;

  //EM algorithm
  do
  {
    V_t_t_s.slice(0) = Q_0; // Q_0 may have been updated or not

    // E-step
    event_time = vecmin(tstart);
#pragma omp parallel                                                                                           \
    private(n_threads_current, i_am, i_points, i_start) \
      firstprivate(U_, u_) default(shared)
      {
        n_threads_current = omp_get_num_threads();
        i_am = omp_get_thread_num();

        for (int t = 1; t < d + 1; t++){ // each thread will have it own
          U_.zeros();
          u_.zeros();

#pragma omp master
          {
            delta_t = I_len[t - 1];
            event_time += delta_t;

            //Rcpp::Rcout << "made "; //TODO: delete me

            // E-step: Filter step
            a_t_less_s.col(t - 1) = F_ *  a_t_t_s.unsafe_col(t - 1);
            V_t_less_s.slice(t - 1) = F_ * V_t_t_s.slice(t - 1) * T_F_ + delta_t * Q;

            // E-step: scoring step: information matrix and scoring vector
            r_set = Rcpp::as<arma::uvec>(risk_sets[t - 1]) - 1;
            n_cols = r_set.size();

            u.zeros();
            U.zeros();

            if(t == d){
              H_diag_inv = arma::vec(n_cols);
              z_dot = arma::mat(n_parems * order_, n_cols);
              z_dot.zeros();
            }

            is_run_parallel = n_cols / 25 > n_threads; // TODO: How to set?
          }

#pragma omp barrier // TODO: is this needed after a omp master?
          if(is_run_parallel || i_am == 0){
            if(is_run_parallel){
              i_points = n_cols / n_threads_current; // size of partition
              i_start = i_am * i_points; // starting array index

              if (i_am == n_threads - 1) // last thread may do more
                i_points = n_cols - i_start;
            } else if(i_am == 0){
              i_start = 0;
              i_points = n_cols;
            }
          }

          if(is_run_parallel || i_am == 0){
            // Get columns to work with
            const arma::uvec i_r_set(r_set.begin() + i_start, i_points, false); // reference the memory
            const arma::vec i_a_t(a_t_less_s.colptr(t - 1), n_parems, false); // reference the memory

            // Compute local result
            unsigned int i = i_start;
            for(auto it = i_r_set.begin(); it != i_r_set.end(); it++){
              const arma::vec x_(_X.colptr(*it), n_parems, false);
              const double i_eta = lower_trunc_exp(arma::dot(i_a_t, x_));

              if(events(*it) && std::abs(tstop(*it) - event_time) < event_eps){
                u_ += x_ * (1.0 - i_eta / (i_eta + 1.0));
              }
              else {
                u_ -= x_ * (i_eta / (i_eta + 1.0));
              }
              U_ += x_ *  (x_.t() * (i_eta / pow(i_eta + 1.0, 2.0))); // I guess this is the fastest http://stackoverflow.com/questions/26766831/armadillo-inplace-plus-significantly-slower-than-normal-plus-operation

              if(t == d){
                H_diag_inv(i) = pow(1.0 + i_eta, 2.0) / i_eta;
                z_dot.rows(0, n_parems - 1).col(i) = x_ *  (i_eta / pow(1.0 + i_eta, 2.0));
                ++i;
              }
            }
          }

          if(is_run_parallel || i_am == 0){
#pragma omp critical
            {
              U.submat(0, 0, n_parems - 1, n_parems - 1) = U.submat(0, 0, n_parems - 1, n_parems - 1) + U_;
              u.head(n_parems) = u.head(n_parems) + u_;
            }
          }


#pragma omp barrier // TODO: is this needed?

#pragma omp master
          {
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

          /* if(t < d + 1 && i_am < 2){ // TODO: Delete
            Rcpp::Rcout << std::setprecision(17) << "It = " << it << " t = " << t << std::endl;
            U.raw_print(Rcpp::Rcout);
            u.raw_print(Rcpp::Rcout);
            V_t_less_s.slice(t - 1).raw_print(Rcpp::Rcout);
            V_t_t_s.slice(t).raw_print(Rcpp::Rcout);
            a_t_t_s.col(t).raw_print(Rcpp::Rcout);
          } */

#pragma omp barrier //TODO is barrier needed after master?
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
