// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins(cpp11)]]
#include <iostream>
#include <thread>
#include <future>

#if defined(USE_OPEN_BLAS) // Used to set the number of threads later
#include "cblas.h"
extern void openblas_set_num_threads(int num_threads);
extern int openblas_get_num_threads();
#define ARMA_USE_OPENBLAS
#define ARMA_DONT_USE_WRAPPER
#else
#define ARMA_USE_BLAS
#endif

// we know these are avialble with all R installations
#define ARMA_USE_LAPACK

#include <RcppArmadillo.h>
#include <armadillo> // has to come after the #define ARMA_DONT_USE_WRAPPER

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
#if defined(NDEBUG)
#undef NDEBUG
#endif

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
inline double relative_norm_change(const arma::vec &prev_est, const arma::vec &new_est){
  return arma::norm(prev_est - new_est, 2) / (arma::norm(prev_est, 2) + 1.0e-10);
}
double (*conv_criteria)(const arma::vec&, const arma::vec&) = relative_norm_change;

// Hepler structure to reference data
// TODO: Figure out what is specific to the EKF
class problem_data{
public:
  // Initalize constants
  const int d;
  const Rcpp::List risk_sets;
  const int n_parems;

  const arma::mat &F_;
  const arma::mat T_F_;

  // Note a copy of data is made below and it is not const for later initalization of pointer to memory (this is not possible with a const pointer)
  arma::mat _X;

  const std::vector<double> I_len;
  const double event_eps; // something small

#if defined(USE_OPEN_BLAS)
  const int n_threads = std::thread::hardware_concurrency();
#endif

  const arma::vec &tstart;
  const arma::vec &tstop;
  const arma::ivec &is_event_in_bin;

  // Declare and maybe intialize non constants
  arma::mat &Q;

  arma::mat a_t_t_s;
  arma::mat a_t_less_s;

  arma::cube V_t_t_s;
  arma::cube V_t_less_s;
  arma::cube B_s;

  arma::cube lag_one_cor;

  problem_data(const Rcpp::NumericMatrix &X, const arma::vec &tstart_,
               const arma::vec &tstop_, const arma::ivec &is_event_in_bin_,
               const arma::colvec &a_0,
               arma::mat &Q_0,
               arma::mat &Q_,
               const Rcpp::List &risk_obj,
               const arma::mat &F__,
               const int n_max = 100, const double eps = 0.001,
               const bool verbose = false, const bool save_all_output = false,
               const int order_ = 1, const bool est_Q_0 = true):
    d(Rcpp::as<int>(risk_obj["d"])),
    risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"])),
    n_parems(a_0.size() / order_),
    F_(F__),
    T_F_(F_.t()),
    I_len(Rcpp::as<std::vector<double> >(risk_obj["I_len"])),
    event_eps(d * std::numeric_limits<double>::epsilon()),
#if defined(USE_OPEN_BLAS)
    n_threads(std::thread::hardware_concurrency()),
#endif
    tstart(tstart_),
    tstop(tstop_),
    is_event_in_bin(is_event_in_bin_),
    Q(Q_)
  {
    // Note a copy of data is made below and it is not const for later initalization of pointer to memory (this is not possible with a const pointer)
    _X = arma::mat(X.begin(), X.nrow(), X.ncol()).t(); // Armadillo use column major ordering https://en.wikipedia.org/wiki/Row-major_order

    a_t_t_s = arma::mat(n_parems * order_, d + 1);
    a_t_less_s = arma::mat(n_parems * order_, d);
    V_t_t_s = arma::cube(n_parems * order_, n_parems * order_, d + 1);
    V_t_less_s = arma::cube(n_parems * order_, n_parems * order_, d);
    B_s = arma::cube(n_parems * order_, n_parems * order_, d);

    a_t_t_s.col(0) = a_0;

    lag_one_cor = arma::cube(n_parems * order_, n_parems * order_, d);
  }
};

// Class with further members for the Extended Kalman Filter
class problem_data_EKF : public problem_data{
public:
  // locks for parallel implementation
  std::mutex m_u;
  std::mutex m_U;

  // Maticies for score and information matrix
  arma::colvec u;
  arma::mat U;

  // Needed for lag one covariance
  arma::mat z_dot;
  arma::vec H_diag_inv;
  arma::mat K_d;

  problem_data_EKF(const Rcpp::NumericMatrix &X, const arma::vec &tstart_,
                 const arma::vec &tstop_, const arma::ivec &is_event_in_bin_,
                 const arma::colvec &a_0,
                 arma::mat &Q_0,
                 arma::mat &Q_,
                 const Rcpp::List &risk_obj,
                 const arma::mat &F__,
                 const int n_max = 100, const double eps = 0.001,
                 const bool verbose = false, const bool save_all_output = false,
                 const int order_ = 1, const bool est_Q_0 = true):
  problem_data(X, tstart_, tstop_, is_event_in_bin_, a_0,
               Q_0, Q_, risk_obj, F__, n_max, eps, verbose, save_all_output,
               order_, est_Q_0)
  {
    u = arma::colvec(n_parems * order_);
    U = arma::mat(n_parems * order_, n_parems * order_);

    const int n_in_last_set = (Rcpp::as<arma::uvec>(risk_sets[d - 1])).size();
    z_dot = arma::mat(n_parems * order_, n_in_last_set);
    H_diag_inv = arma::vec(n_in_last_set);
  }
};

// Abstact solver class
class Solver {
public:
  virtual void solve() = 0;
};












// Class to handle parallel computation in extended Kalman filter
class EKF_helper{
  using uvec_iter = arma::uvec::const_iterator;

  // handy class to ensure that all ressources are cleaned up once scope of the
  // class is expired
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
  // This class is abstact as the method do_computation will differ between
  // the models
  class filter_worker{
  protected:
    virtual void do_comps(const uvec_iter it, int &i,
                          const arma::vec &i_a_t, const bool &compute_z_and_H,
                          const double &event_time, const int &bin_number)
    {}; // abstact method to be implemented

    bool is_first_call;
    problem_data_EKF &dat;

    // local variables to compute temporary result
    arma::colvec u_;
    arma::mat U_;

  public:
    filter_worker(problem_data_EKF &p_data):
    is_first_call(true), dat(p_data)
    {}

    bool operator()(uvec_iter first, uvec_iter last,
                    const arma::vec &i_a_t, bool compute_z_and_H,
                    double event_time, int i_start, const int bin_number){
      // potentially intialize variables and set entries to zeroes in any case
      if(is_first_call){
        u_ = arma::vec(dat.n_parems);
        U_ = arma::mat(dat.n_parems, dat.n_parems);
        is_first_call = false;
      }
      u_.zeros();
      U_.zeros();

      /*{ //TODO: delete
        std::lock_guard<std::mutex> lk(dat.m_u);
        Rcpp::Rcout << "I am " << std::this_thread::get_id() << " " <<
          std::distance(first, last) << " " <<
            compute_z_and_H << " " <<  event_time << " " << bin_number << std::endl;
        i_a_t.print();
      }*/

      // compute local results
      int i = i_start;
      for(uvec_iter it = first; it != last; it++){
        do_comps(it, i, i_a_t, compute_z_and_H, event_time, bin_number);
      }

      // Update shared variable
      {
        std::lock_guard<std::mutex> lk(dat.m_U);
        dat.U.submat(0, 0, dat.n_parems - 1, dat.n_parems - 1) +=  U_;
      }

      {
        std::lock_guard<std::mutex> lk(dat.m_u);
        /*Rcpp::Rcout << std::this_thread::get_id() << " got this to add in bin " << bin_number << std::endl;
        u_.print();

        Rcpp::Rcout << "Here is u_ before" << std::endl;

        dat.u.print();*/

        dat.u.head(dat.n_parems) += u_;

        /*Rcpp::Rcout << "Here is the new u_" << std::endl;
        dat.u.print();

        Rcpp::Rcout << std::endl;*/
      }

      return true;
    }
  };

  class filter_worker_logit : public filter_worker {
  private:
    void do_comps(const uvec_iter it, int &i,
                  const arma::vec &i_a_t, const bool &compute_z_and_H,
                  const double &event_time, const int &bin_number){
      const arma::vec x_(dat._X.colptr(*it), dat.n_parems, false);
      const double i_eta = lower_trunc_exp(arma::dot(i_a_t, x_));

      if(dat.is_event_in_bin(*it) == bin_number){
        u_ += x_ * (1.0 - i_eta / (i_eta + 1.0));
      }
      else {
        u_ -= x_ * (i_eta / (i_eta + 1.0));
      }
      U_ += x_ *  (x_.t() * (i_eta / pow(i_eta + 1.0, 2.0))); // I guess this is the fastest http://stackoverflow.com/questions/26766831/armadillo-inplace-plus-significantly-slower-than-normal-plus-operation

      /*{ //TODO: delete
        std::lock_guard<std::mutex> lk(dat.m_u);
        Rcpp::Rcout << "I am " << std::this_thread::get_id() << " boh!" << std::endl;
        (x_ *  (x_.t() * (i_eta / pow(i_eta + 1.0, 2.0)))).print();
      }*/

      if(compute_z_and_H){
        dat.H_diag_inv(i) = pow(1.0 + i_eta, 2.0) / i_eta;
        dat.z_dot.rows(0, dat.n_parems - 1).col(i) = x_ *  (i_eta / pow(1.0 + i_eta, 2.0));
        ++i;
      }

      /*std::stringstream ss; TODO: delete
      ss << std::this_thread::get_id() << "\t" << &u_ << "\t" << &dat.m_u << "\t" << &dat.U << "\n";
      Rcpp::Rcout << ss.str();*/
    }

  public:
    filter_worker_logit(problem_data_EKF &p_data):
      filter_worker(p_data)
    {}
  };

  unsigned long const hardware_threads;
  problem_data_EKF &p_data;
  std::vector<std::shared_ptr<filter_worker> > workers;

public:
  EKF_helper(problem_data_EKF &p_data_):
  hardware_threads(std::thread::hardware_concurrency()),
  p_data(p_data_), workers()
  {
    // create workers
    for(int i = 0; i < hardware_threads; i++){
      std::shared_ptr<filter_worker> new_p(new filter_worker_logit(p_data));
      workers.push_back(std::move(new_p));
    }
  }

  void parallel_filter_step(uvec_iter first, uvec_iter last,
                            const arma::vec &i_a_t,
                            const bool compute_H_and_z,
                            double event_time, const int bin_number){
    // Set referenced objects entries to zero
    p_data.U.zeros();
    p_data.u.zeros();
    if(compute_H_and_z){
      p_data.z_dot.zeros();
      p_data.H_diag_inv.zeros();
    }

    // Compute the number of threads to create
    unsigned long const length = std::distance(first, last);

    unsigned long const min_per_thread = 50; // TODO: figure out how to set?
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
      std::packaged_task<bool(
        uvec_iter, uvec_iter, const arma::vec&, bool, double, int, const int)> task(
            [=](uvec_iter block_start, uvec_iter block_end, const arma::vec& i_a_t, bool compute_H_and_z,
               double event_time, int i_start, const int bin_number){
              return (*it->get())(block_start, block_end, i_a_t, compute_H_and_z,
                      event_time, i_start, bin_number);
            });
      futures[i] = task.get_future();
      threads[i] = std::thread(std::move(task), block_start, block_end, i_a_t, compute_H_and_z,
                               event_time, i_start, bin_number);

      i_start += block_size;
      block_start = block_end;
    }
    (*(it->get()))(block_start, last, i_a_t, compute_H_and_z, event_time, i_start, bin_number); // compute last enteries on this thread

    for(unsigned long i = 0; i < num_threads - 1; ++i)
    {
      futures[i].get(); // will throw if any of the threads did
    }
  }
};

class EKF_solver : public Solver{
  problem_data_EKF &p_dat;
  EKF_helper filter_helper;
  const double min_time;


public:
  EKF_solver(problem_data_EKF &p_):
    p_dat(p_), filter_helper(EKF_helper(p_)),
    min_time(vecmin(p_dat.tstart))
  {}

  void solve(){
    double event_time = min_time;
    for (int t = 1; t < p_dat.d + 1; t++){
      double delta_t = p_dat.I_len[t - 1];
      event_time += delta_t;

      // E-step: Filter step
      p_dat.a_t_less_s.col(t - 1) = p_dat.F_ *  p_dat.a_t_t_s.unsafe_col(t - 1);
      p_dat.V_t_less_s.slice(t - 1) = p_dat.F_ * p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ + delta_t * p_dat.Q;

      // E-step: scoring step: information matrix and scoring vector
      arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;
      const arma::vec i_a_t(p_dat.a_t_less_s.colptr(t - 1), p_dat.n_parems, false);

#if defined(USE_OPEN_BLAS)
      openblas_set_num_threads(1);
      //Rcpp::Rcout << "n thread before = " << openblas_get_num_threads() << std::endl;
#endif

      filter_helper.parallel_filter_step(r_set.begin(), r_set.end(), i_a_t, t == p_dat.d, event_time, t - 1);

#ifdef USE_OPEN_BLAS
      openblas_set_num_threads(p_dat.n_threads);
      //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif

      // E-step: scoring step: update values
      arma::mat V_t_less_s_inv = arma::inv_sympd(p_dat.V_t_less_s.slice(t - 1));
      p_dat.V_t_t_s.slice(t) = arma::inv_sympd(V_t_less_s_inv + p_dat.U);
      p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.unsafe_col(t - 1) + p_dat.V_t_t_s.slice(t) * p_dat.u;
      p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * V_t_less_s_inv;

      if(t == p_dat.d){
        p_dat.K_d = p_dat.V_t_less_s.slice(t - 1) * inv(arma::eye<arma::mat>(size(p_dat.U)) + p_dat.U *
          p_dat.V_t_less_s.slice(t - 1)) * p_dat.z_dot * diagmat(p_dat.H_diag_inv);
        // Parenthesis is key here to avoid making a n x n matrix for large n
        p_dat.K_d = (p_dat.F_ * p_dat.V_t_less_s.slice(t - 1) * p_dat.z_dot * diagmat(p_dat.H_diag_inv) * p_dat.z_dot.t()) * p_dat.K_d;
        p_dat.K_d = p_dat.F_ * p_dat.V_t_less_s.slice(t - 1) * p_dat.z_dot * diagmat(p_dat.H_diag_inv) -  p_dat.K_d;

        p_dat.lag_one_cor.slice(t - 1) = (arma::eye<arma::mat>(size(p_dat.U)) - p_dat.K_d * p_dat.z_dot.t()) * p_dat.F_ * p_dat.V_t_t_s.slice(t - 1);
      }
    }
  }
};






class UKF_solver : public Solver{

  problem_data &p_dat;
  const arma::uword m;
  const double k;
  const double w_0;
  const double w_i;
  const double min_time;
  const double sqrt_m_k;
  arma::mat sigma_points;

  inline void compute_sigma_points(const arma::vec &a_t,
                                   arma::mat &s_points,
                                   const arma::mat &P_x_x){
    const arma::mat cholesky_decomp = arma::chol(P_x_x, "upper").t(); // TODO: cholesky_decomp * cholesky_decomp.t() = inital mat. should ensure that . See http://arma.sourceforge.net/docs.html#chol

    s_points.col(0) = a_t;
    for(int i = 1; i < s_points.n_cols; ++i)
      if(i % 2 == 0)
        s_points.col(i) = a_t + sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2); else
          s_points.col(i) = a_t - sqrt_m_k * cholesky_decomp.unsafe_col((i - 1) / 2);
  }

public:
  UKF_solver(problem_data &p_, Rcpp::Nullable<Rcpp::NumericVector> &k_):
    p_dat(p_),
    m(p_.a_t_t_s.n_rows),
    k(!k_.isNull() ? Rcpp::as< Rcpp::NumericVector >(k_)[0] : 3.0 - m),
    w_0(k / (m + k)),
    w_i(1 / (2 * (m + k))),
    min_time(vecmin(p_dat.tstart)),
    sqrt_m_k(std::sqrt(m + k)),
    sigma_points(arma::mat(m, 2 * m + 1))
    {}

  void solve(){
#ifdef USE_OPEN_BLAS //TODO: Move somewhere else?
    openblas_set_num_threads(p_dat.n_threads);
    //Rcpp::Rcout << "n thread after = " << openblas_get_num_threads() << std::endl;
#endif

    double event_time = min_time;

    // Need to compute the first set of sigma points
    compute_sigma_points(p_dat.a_t_t_s.unsafe_col(0),
                         sigma_points, p_dat.V_t_t_s.slice(0));

    for (int t = 1; t < p_dat.d + 1; t++){
      double delta_t = p_dat.I_len[t - 1];
      event_time += delta_t;

      // E-step: Filter step
      //   Updates a_t_less_s and V_t_less_s
      //   Requires T(sigma point) + a_t_less_s computed before V_t_less_s
      //     NB have to compute sigma points the first time arround
      //     Requires for-loop with 2m + 1 itertions

      // First we compute the mean
      p_dat.a_t_less_s.col(t - 1) = w_0 * sigma_points.unsafe_col(0) +
        w_i * arma::sum(sigma_points.cols(1, sigma_points.n_cols - 1), 1);

      // Then the variance
      p_dat.V_t_less_s.slice(t - 1) = p_dat.Q; // weigths sum to one // TODO: Include or not?
      for(int i = 0; i < sigma_points.n_cols; ++i){
        const double &w = i == 0 ? w_0 : w_i;

        p_dat.V_t_less_s.slice(t - 1) +=
          (w * (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
            (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1)).t();
      }

      // Recompute sigma points
      //TODO: should or should these be re-computed?
      //compute_sigma_points(p_dat.a_t_less_s.unsafe_col(t - 1),
      //                     sigma_points, p_dat.V_t_less_s.slice(t - 1));

      // E-step: correction-step
      //   Compute a_t_t_s and v_t_t_s
      //   Update y_bar, P_a_v, P_v_v
      //     Can be done in 2m + 1 iterations
      //   Update a_t_t_s
      //   Requires comptation of the (2m + 1) x | risk_set |  matrix
      //   Need to update sigma points (TODO: unless it is the last iteration?)
      //     Thus a Cholesky decomposition V_t_less_s

      // First, we compute the mean of the outcome
      arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;

      arma::mat &&tmp = (sigma_points.t() * p_dat._X.cols(r_set)).t(); // we transpose due to the column-major
      arma::mat Z_t = tmp.for_each(
        [](arma::mat::elem_type &val) {
          double &&tmp = lower_trunc_exp(val);
          val = tmp / (1 + tmp);
          });

      // Compute y_bar, P_a_v and P_v_v
      const arma::vec y_bar = w_0 * Z_t.unsafe_col(0) +
        w_i * arma::sum(Z_t.cols(1, Z_t.n_cols - 1), 1);

      arma::mat P_a_v = (w_0 * (sigma_points.unsafe_col(0) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
        (Z_t.unsafe_col(0) - y_bar).t();

      arma::mat P_v_v = (w_0 * (Z_t.unsafe_col(0) - y_bar)) * (Z_t.unsafe_col(0) - y_bar).t() +
        //arma::diagmat(y_bar % (1 - y_bar));
        arma::diagmat(w_0 * Z_t.unsafe_col(0) % (1 - Z_t.unsafe_col(0))); // % is element-wise product. TODO: I think we have to time with w_0 here too? same applies for w_i below

      for(int i = 1; i < sigma_points.n_cols; ++i){
        P_a_v += (w_i * (sigma_points.unsafe_col(i) - p_dat.a_t_less_s.unsafe_col(t - 1))) *
          (Z_t.unsafe_col(i) - y_bar).t();
        P_v_v += (w_i * (Z_t.unsafe_col(i) - y_bar)) * (Z_t.unsafe_col(i) - y_bar).t() +
          arma::diagmat(w_i * Z_t.unsafe_col(i) % (1 - Z_t.unsafe_col(i)));
      }

      // Compute new estimates
      P_v_v = arma::inv(P_v_v); // NB: Note that we invert the matrix here so P_v_v is inv(P_v_v)

      p_dat.a_t_t_s.col(t) = p_dat.a_t_less_s.unsafe_col(t - 1) +
        P_a_v * (P_v_v * (
            (p_dat.is_event_in_bin(r_set) == t -1) - y_bar));

      p_dat.V_t_t_s.slice(t) = p_dat.V_t_less_s.slice(t - 1) - P_a_v * P_v_v * P_a_v.t();

      p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * arma::inv_sympd(p_dat.V_t_less_s.slice(t - 1));

      // Update sigma pooints for next iteration
      compute_sigma_points(p_dat.a_t_t_s.unsafe_col(t),
                           sigma_points, p_dat.V_t_t_s.slice(t));
    }
  }
};




//' @export
// [[Rcpp::export]]
Rcpp::List ddhazard_fit_cpp_prelim(const Rcpp::NumericMatrix &X, const arma::vec &tstart,
                                   const arma::vec &tstop,
                                   const arma::colvec &a_0,
                                   arma::mat Q_0, // by value copy. This  is key cuz we will change it if est_Q_0 = T
                                   arma::mat Q, // similary this is a copy
                                   const Rcpp::List &risk_obj,
                                   const arma::mat &F_,
                                   const int n_max = 100, const double eps = 0.001,
                                   const bool verbose = false, const bool save_all_output = false,
                                   const int order_ = 1, const bool est_Q_0 = true,
                                   const std::string method = "EKF",
                                   Rcpp::Nullable<Rcpp::NumericVector> k = R_NilValue // see this link for nullable example http://blogs.candoerz.com/question/164706/rcpp-function-for-adding-elements-of-a-vector.aspx
                                     ){
  // Declare and maybe intialize non constants
  double delta_t, test_max_diff;
  const double Q_warn_eps = sqrt(std::numeric_limits<double>::epsilon());

  Rcpp::NumericVector conv_values;

  arma::colvec a_prev(a_0.begin(), a_0.size());
  Rcpp::List all_output; // only need if save_all_output = true

  unsigned int it = 0;

  // M-stp pointers for convenience
  // see http://stackoverflow.com/questions/35983814/access-column-of-matrix-without-creating-copy-of-data-using-armadillo-library
  // the use of unsafe_col is key
  arma::mat *B, *V_less, *V;
  arma::vec a_less, a;

  const arma::ivec is_event_in_bin = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);

  // Intialize the solver for the E-step
  problem_data *p_data;
  Solver  *solver;

  if(method == "EKF"){
    p_data = new problem_data_EKF(
      X, tstart, tstop, is_event_in_bin,
      a_0, Q_0, Q,
      risk_obj, F_,
      n_max, eps, verbose, save_all_output,
      order_, est_Q_0);
    solver = new EKF_solver(static_cast<problem_data_EKF &>(*p_data));
  } else if (method == "UKF"){
    p_data = new problem_data(
      X, tstart, tstop, is_event_in_bin,
      a_0, Q_0, Q,
      risk_obj, F_,
      n_max, eps, verbose, save_all_output,
      order_, est_Q_0);
    solver = new UKF_solver(*p_data, k);
  } else
    Rcpp::stop("method '" + method  +"'is not implemented");

  /*Rcpp::Rcout << "X" << std::endl; TODO: Clean up or move
  p_data._X.print();
  Rcpp::Rcout << "a_t_less_s" << std::endl;
  p_data.a_t_less_s.print();
  Rcpp::Rcout << "a_t_t_s" << std::endl;
  p_data.a_t_t_s.print();
  Rcpp::Rcout << "B_s" << std::endl;
  p_data.B_s.print();
  Rcpp::Rcout << "d " << p_data.d << std::endl;
  Rcpp::Rcout << "event_eps " << p_data.event_eps << std::endl;
  Rcpp::Rcout << "is_event_in_bin" << std::endl;
  p_data.is_event_in_bin.print();
  Rcpp::Rcout << "F_" << std::endl;
  p_data.F_.print();
  Rcpp::Rcout << "H_diag_inv" << std::endl;
  p_data.H_diag_inv.print();
  Rcpp::Rcout << "I_len size" << p_data.I_len.size() << std::endl;
  Rcpp::Rcout << "K_d" << std::endl;
  p_data.K_d.print();
  Rcpp::Rcout << "lag_one_cor" << std::endl;
  p_data.lag_one_cor.print();
  Rcpp::Rcout << "risk_sets size" << p_data.risk_sets.size() <<  std::endl;
  Rcpp::Rcout << "T_F_" << std::endl;
  p_data.T_F_.print();
  Rcpp::Rcout << "tstart" << std::endl;
  p_data.tstart.print();
  Rcpp::Rcout << "tstop" << std::endl;
  p_data.tstop.print();
  Rcpp::Rcout << "U" << std::endl;
  p_data.U.print();
  Rcpp::Rcout << "u" << std::endl;
  p_data.u.print();
  Rcpp::Rcout << "V_t_less_s" << std::endl;
  p_data.V_t_less_s.print();
  Rcpp::Rcout << "V_t_t_s" << std::endl;
  p_data.V_t_t_s.print();
  Rcpp::Rcout << "z_dot" << std::endl;
  p_data.z_dot.print();*/

  do
  {
    if((it + 1) % 25 == 0)
      Rcpp::checkUserInterrupt(); // this is expensive - you do not want to check too often

    p_data->V_t_t_s.slice(0) = Q_0; // Q_0 may have been updated or not

    // E-step
    solver->solve();

    // E-step: smoothing
    for (int t = p_data->d - 1; t > -1; t--){
      // we need to compute the correlation matrix first
      if(t > 0){
          p_data->lag_one_cor.slice(t - 1) = p_data->V_t_t_s.slice(t) * p_data->B_s.slice(t - 1).t() +
          p_data->B_s.slice(t) * (
          p_data->lag_one_cor.slice(t) - F_ * p_data->V_t_t_s.slice(t)) * p_data->B_s.slice(t - 1).t();
      }

      p_data->a_t_t_s.col(t) = p_data->a_t_t_s.unsafe_col(t) + p_data->B_s.slice(t) *
        (p_data->a_t_t_s.unsafe_col(t + 1) - p_data->a_t_less_s.unsafe_col(t));
      p_data->V_t_t_s.slice(t) = p_data->V_t_t_s.slice(t) + p_data->B_s.slice(t) *
        (p_data->V_t_t_s.slice(t + 1) - p_data->V_t_less_s.slice(t)) * p_data->B_s.slice(t).t();
  }

    // M-step
    if(est_Q_0){
      Q_0 = p_data->V_t_t_s.slice(0);
    }
    Q.zeros();
    for (int t = 1; t < p_data->d + 1; t++){
      delta_t = p_data->I_len[t - 1];

      B = &p_data->B_s.slice(t - 1);
      V_less = &p_data->V_t_t_s.slice(t - 1);
      V = &p_data->V_t_t_s.slice(t);
      a_less = p_data->a_t_t_s.unsafe_col(t - 1);
      a = p_data->a_t_t_s.unsafe_col(t);

      Q += ((a - F_ * a_less) * (a - F_ * a_less).t() + *V
              - F_ * *B * *V
              - (F_ * *B * *V).t()
              + F_ * *V_less * p_data->T_F_) / delta_t;
    }
    Q /= p_data->d;

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
      arma::mat tmp_Q = Q.submat(0, 0, p_data->n_parems - 1, p_data->n_parems - 1);
      Q.zeros();
      Q.submat(0, 0, p_data->n_parems - 1, p_data->n_parems - 1) = tmp_Q;
    }

    conv_values.push_back(conv_criteria(a_prev, p_data->a_t_t_s.unsafe_col(0)));

    //if(save_all_output) // TODO: make similar save all output function?

    if(*(conv_values.end() -1) < eps){
      break;
    }

    if(verbose && it % 5 == 0){
      Rcpp::Rcout << "Iteration " << it + 1 << " ended with conv criteria " << *(conv_values.end() -1) << std::endl;
    }

    a_prev = p_data->a_t_t_s.col(0);
}while(++it < n_max);

  /* if(save_all_output){
  all_output
  }else */

  if(it == n_max)
    throw std::runtime_error("EM algorithm did not converge within the n_max number of iterations");

  return(Rcpp::List::create(Rcpp::Named("V_t_d_s") = Rcpp::wrap(p_data->V_t_t_s),
                            Rcpp::Named("a_t_d_s") = Rcpp::wrap(p_data->a_t_t_s.t()),
                            Rcpp::Named("B_s") = Rcpp::wrap(p_data->B_s),
                            Rcpp::Named("lag_one_cor") = Rcpp::wrap(p_data->lag_one_cor),

                            Rcpp::Named("n_iter") = it + 1,
                            Rcpp::Named("conv_values") = conv_values,
                            Rcpp::Named("Q") = Rcpp::wrap(Q),
                            Rcpp::Named("Q_0") = Rcpp::wrap(Q_0)));
  }


/*** R
library(survival); library(testthat); library(dynamichazard); source("../R/test_utils.R")

set.seed(2972)
sims <- test_sim_func_logit(n_series = 10^4, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = .1, x_mean = -.4, re_draw = T)
sims$res <- as.data.frame(sims$res)

design_mat <- get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3 - 1, sims$res)
rist_sets <- get_risk_obj(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

log_file = file("debug.log")
sink(log_file)

arg_list <- list(
  X = design_mat$X,
  tstart = design_mat$Y[, 1],  tstop = design_mat$Y[, 2],
  a_0 = rep(0, ncol(design_mat$X)),
  Q_0 = diag(10, ncol(design_mat$X)), # something large
  Q = diag(1, ncol(design_mat$X)), # something large
  F_ = diag(1, ncol(design_mat$X)), # first order random walk
  risk_obj = rist_sets,
  eps = 10^-4, n_max = 10^4,
  order_ = 1,
  est_Q_0 = F
)

res <- do.call(ddhazard_fit, arg_list)

tryCatch({
  res_new <- do.call(ddhazard_fit_cpp_prelim, arg_list)
}, finally = {
  sink()
  close(log_file)
})


test_that("Testing old versus new implementation", {
  expect_equal(res$a_t_d_s, res_new$a_t_d_s)
  expect_equal(res$V_t_d_s, res_new$V_t_d_s)
  expect_equal(res$B_s, res_new$B_s)
})

# Test UKF
sims <- test_sim_func_logit(n_series = 10^3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = .1, x_mean = -.4, re_draw = T)
sims$res <- as.data.frame(sims$res)

log_file = file("debug.log", open = "a")
sink(log_file)

design_mat <- get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3 - 1, sims$res)
sum(design_mat$Y[, 3] & design_mat$Y[, 2] <= 10)
rist_sets <- get_risk_obj(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

design_mat$Y[, 2] <- rist_sets$stop_new
design_mat$Y[, 3] <- rist_sets$new_events_flags

arg_list <- list(
  X = design_mat$X,
  tstart = design_mat$Y[, 1],  tstop = design_mat$Y[, 2], events = design_mat$Y[, 3],
  a_0 = c(1, 1, 1),
  Q_0 = diag(1.0e1, ncol(design_mat$X)), # something large
  Q = diag(1.0e-0, ncol(design_mat$X)), # something large
  F_ = diag(1, ncol(design_mat$X)), # first order random walk
  risk_obj = rist_sets,
  eps = 10^-2, n_max = 10^4,
  order_ = 1,
  est_Q_0 = F,
  verbose = T,
  method = "UKF"
)

tryCatch({
  cat("UKF runtime\n")
  print(system.time(res_UKF <- do.call(ddhazard_fit_cpp_prelim, arg_list)))
  arg_list$method <- "EKF"
  cat("EKF runtime\n")
  print(system.time(res_EKF <- do.call(ddhazard_fit_cpp_prelim, arg_list)))

  cat("MSE for UKF\n")
  print(mean((res_UKF$a_t_d_s - sims$betas)^2))
  cat("MSE for EKF\n")
  print(mean((res_EKF$a_t_d_s - sims$betas)^2))
}, finally = {
  sink()
  close(log_file)
})
*/
