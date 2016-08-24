// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <omp.h>
#include <iostream>
#include <armadillo>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

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
#define ARMA_DEFAULT_OSTREAM Rcout
#endif

// R now defines NDEBUG which suppresses a number of useful
// Armadillo tests Users can still defined it later, and/or
// define ARMA_NO_DEBUG
#if defined(NDEBUG)
#undef NDEBUG
#endif

// Maybe look at this one day https://github.com/Headtalk/armadillo-ios/blob/master/armadillo-4.200.0/include/armadillo_bits/fn_trunc_exp.hpp
const double lower_trunc_exp_log_thres = sqrt(log(std::numeric_limits<double>::max())) - 1.1;
const double lower_trunc_exp_exp_thres = exp(lower_trunc_exp_log_thres);
void in_place_lower_trunc_exp(arma::colvec & result)
{
  double * x;

  for(unsigned i = 0; i < result.size(); i++){
    x = &result(i);
    if(*x >= lower_trunc_exp_log_thres )
    {
      *x =  lower_trunc_exp_exp_thres;
    }
    else
    {
      *x = std::exp(*x);
    }
  }
}

// from http://gallery.rcpp.org/articles/vector-minimum/
template<typename T>
int vecmin(const T x){
  // Rcpp supports STL-style iterators
  auto it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}

//' Function to perform one iteration of the Kalman filter in parallel
//' @export
// [[Rcpp::export]]
List gen_kalman_filter_cpp(const arma::colvec &a_0, const arma::mat &Q_0, const arma::mat &Q, const arma::mat &F_,
                           const List &risk_sets, const NumericVector &I_len,
                           const int d, const NumericMatrix &X,
                           const IntegerVector &start, const arma::ivec &stop, const arma::ivec &events, // Armadillo have no boolean vector - do not use uvec fasting fails http://stackoverflow.com/questions/32309915/rcpparmadillo-c-create-bool-vector
                           const int & order_){
  // see http://stackoverflow.com/questions/18866130/passing-large-matrices-to-rcpparmadillo-function-without-creating-copy-advanced
  int event_time, delta_t;

  const arma::mat T_F_ = F_.t();
  const arma::mat _X = arma::mat(X.begin(), X.nrow(), X.ncol()).t(); // Armadillo use column major ordering https://en.wikipedia.org/wiki/Row-major_order
  const int n_parems = a_0.size() / order_;

  // Define matrices for later
  arma::mat a_t_t_s = arma::mat(n_parems * order_, d + 1);
  arma::mat a_t_less_s = arma::mat(n_parems * order_, d);

  // Last input is the slice argument
  cube V_t_t_s = cube(n_parems * order_, n_parems * order_, d + 1);
  cube V_t_less_s = cube(n_parems * order_, n_parems * order_, d);
  cube B_s = cube(n_parems * order_, n_parems * order_, d);

  a_t_t_s.col(0) = a_0;
  V_t_t_s.slice(0) = Q_0;

  arma::colvec u = arma::colvec(n_parems * order_);
  arma::mat U = arma::mat(n_parems * order_, n_parems * order_);

  uvec r_set;
  arma::colvec exp_eta;
  arma::mat V_t_less_s_inv;
  double exp_eta_it;

  // Needed for lag one covariance
  arma::mat z_dot;
  vec H_diag_inv;
  arma::mat K_d;
  cube lag_one_cor = cube(n_parems * order_, n_parems * order_, d);

  // Parallel settings and variables
  const int n_threads = omp_get_num_procs() - 1;
  omp_set_num_threads(n_threads);
  int i_am, i_points, i_start, n_cols, n_threads_current;
  arma::mat i_x_;
  uvec i_r_set;
  arma::ivec i_stop;
  arma::ivec i_events;

  unsigned int i;
  arma::colvec u_ = arma::colvec(n_parems);
  arma::mat U_ = arma::mat(n_parems, n_parems);
  bool is_run_parallel;

  // Start loop
  event_time = vecmin(start);
  // See http://stackoverflow.com/questions/31321071/openmp-nested-for-loop-becomes-faster-when-having-parallel-before-outer-loop
#pragma omp parallel                                                                                           \
  private(n_threads_current, exp_eta_it, i_am, i_points, i_start, i_x_, exp_eta, i, i_r_set, i_stop, i_events) \
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

  // Filter step
  a_t_less_s.col(t - 1) = F_ *  a_t_t_s.col(t - 1);
  V_t_less_s.slice(t - 1) = F_ * V_t_t_s.slice(t - 1) * T_F_ + delta_t * Q;

  // scoring step: information matrix and scoring vector
  r_set = as<uvec>(risk_sets[t - 1]) - 1;
  n_cols = r_set.size();

  u.zeros();
  U.zeros();

  if(t == d){
    H_diag_inv = vec(n_cols);
    z_dot = arma::mat(n_parems * order_, n_cols);
    z_dot.zeros();
  }

  is_run_parallel = n_cols / 3 > n_threads; // TODO: How to set?
}

#pragma omp barrier // TODO: is this needed after a omp master?

if(is_run_parallel){
  i_points = n_cols / n_threads_current; // size of partition
  i_start = i_am * i_points; // starting array index

  if (i_am == n_threads - 1) // last thread may do more
    i_points = n_cols - i_start;

  // Get columns to work with
  i_r_set = r_set(span(i_start, i_start + i_points - 1));
  i_x_  = _X.cols(i_r_set); // This is not reference copy but value http://stackoverflow.com/questions/18859328/fastest-way-to-refer-to-vector-in-armadillo-library
  exp_eta =  i_x_.t() * a_t_less_s.col(t - 1).head(n_parems);
  in_place_lower_trunc_exp(exp_eta);

  i_events = events(i_r_set);
  i_stop = stop(i_r_set);

  if(t == d){
    H_diag_inv(span(i_start, i_start + i_points - 1)) = pow(1.0 + exp_eta, 2) / exp_eta;
    z_dot.rows(0, n_parems - 1).cols(i_start, i_start + i_points - 1) =
      i_x_ *  diagmat(exp_eta / pow(1.0 + exp_eta, 2));
  }
}
else if (i_am == 0){
  i_x_  = _X.cols(r_set); // This is not reference copy but value http://stackoverflow.com/questions/18859328/fastest-way-to-refer-to-vector-in-armadillo-library
  exp_eta =  i_x_.t() * a_t_less_s.col(t - 1).head(n_parems);
  in_place_lower_trunc_exp(exp_eta);

  i_events = events(r_set);
  i_stop = stop(r_set);

  if(t == d){
    H_diag_inv = pow(1.0 + exp_eta, 2)/exp_eta;
    z_dot.rows(0, n_parems - 1) = i_x_ *  diagmat(exp_eta / pow(1.0 + exp_eta, 2));
  }
}

if(is_run_parallel || i_am ==0){
  // Compute local result
  for(i = 0; i < exp_eta.size(); i++){
    exp_eta_it = exp_eta(i);
    if(i_events(i) && i_stop(i) == event_time){
      u_ = u_ + i_x_.col(i) * (1.0 - exp_eta_it / (exp_eta_it + 1.0));
    }
    else {
      u_ = u_ - i_x_.col(i) * exp_eta_it / (exp_eta_it + 1.0);
    }
    U_ = U_ + i_x_.col(i) * i_x_.col(i).t() * exp_eta_it / pow(exp_eta_it + 1.0, 2.0);
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
  // scoring step: update values
  V_t_less_s_inv = inv_sympd(V_t_less_s.slice(t - 1));
  V_t_t_s.slice(t) = inv_sympd(V_t_less_s_inv + U);
  a_t_t_s.col(t) = a_t_less_s.col(t - 1) + V_t_t_s.slice(t) * u;
  B_s.slice(t - 1) = V_t_t_s.slice(t - 1) * T_F_ * V_t_less_s_inv;

  if(t == d){
    K_d = V_t_less_s.slice(t - 1) * inv(eye<arma::mat>(size(U)) + U * V_t_less_s.slice(t - 1)) * z_dot * diagmat(H_diag_inv);
    // Parenthesis is key here to avoid making a n x n matrix for large n
    K_d = (F_ * V_t_less_s.slice(t - 1) * z_dot * diagmat(H_diag_inv) * z_dot.t()) * K_d;
    K_d = F_ * V_t_less_s.slice(t - 1) * z_dot * diagmat(H_diag_inv) -  K_d;

    lag_one_cor.slice(t - 1) = (eye<arma::mat>(size(U)) - K_d * z_dot.t()) * F_ * V_t_t_s.slice(t - 1);
  }
}

#pragma omp barrier //TODO is barrier needed after master?
      }
    }

  // smoothing
  for (int t = d - 1; t > -1; t--){
    // we need to compute the correlation matrix first
    if(t > 0){
      lag_one_cor.slice(t - 1) = V_t_t_s.slice(t) * B_s.slice(t - 1).t() +  B_s.slice(t) * (
        lag_one_cor.slice(t) - F_ * V_t_t_s.slice(t)) * B_s.slice(t - 1).t();
    }

    a_t_t_s.col(t) = a_t_t_s.col(t) + B_s.slice(t) *
      (a_t_t_s.col(t + 1) - a_t_less_s.col(t));
    V_t_t_s.slice(t) = V_t_t_s.slice(t) + B_s.slice(t) *
      (V_t_t_s.slice(t + 1) - V_t_less_s.slice(t)) * B_s.slice(t).t();
  }

  return(List::create(Named("V_t_d_s") = wrap(V_t_t_s),
                      Named("a_t_d_s") = wrap(a_t_t_s.t()),
                      Named("B_s") = wrap(B_s),
                      Named("lag_one_cor") = wrap(lag_one_cor)));
}
