/*
 Similar function glm.fit in r-source/src/library/stats/R/glm.R
 It only returns the coefficient vector
*/

#include "arma_n_rcpp.h"
#include "thread_pool.h"
#include "arma_BLAS_LAPACK.h"

template<typename family>
class parallelglm_class {
  using uword = arma::uword;
  static const uword block_size = 1000;

  /* data holder class */
  struct data_holder {
    arma::mat &X;
    arma::vec *beta;
    const arma::vec &Ys;
    const arma::vec &weights;
    const arma::vec &offsets;

    const uword max_threads, p, n;

    arma::mat XtopWX;
    arma::vec XtopWz;

    std::mutex m_XtopWX;
    std::mutex m_XtopWz;

    data_holder(
      arma::mat &X, const arma::vec &Ys, const arma::vec &weights,
      const arma::vec &offsets, const uword max_threads,
      uword p, uword n):
      X(X), Ys(Ys), weights(weights), offsets(offsets), max_threads(max_threads),
      p(p), n(n)
    {}
  };

  /* worker class for multithreading*/
  class worker{
    const uword i_start, i_end;
    data_holder &data;
    const bool first_it;

  public:
    worker(uword i_start, uword i_end, data_holder &data, bool first_it):
    i_start(i_start), i_end(i_end), data(data), first_it(first_it) {}

    void operator()(){

      arma::span my_span(i_start, i_end);
      arma::mat my_X(data.X.begin() + i_start * data.p,
                     data.p, i_end - i_start + 1,
                     false /* dont take copy */);
      arma::vec etas;
      if(first_it){
        etas = arma::vec(i_end - i_start + 1);
        double *eta = etas.begin();
        const double *y = data.Ys.begin() + i_start;
        const double *weight = data.weights.begin() + i_start;
        for(uword i = 0; i < etas.n_elem; ++i, ++eta, ++weight, ++y)
          *eta = family::initialize(*y, *weight);

      } else
        etas = (data.beta->t() * my_X).t() + data.offsets(my_span);

      arma::mat my_XtopWX(data.p, data.p, arma::fill::zeros);
      arma::vec my_XtopWz = arma::vec(data.p, arma::fill::zeros);

      const double *eta = etas.begin();
      const double *weight = data.weights.begin() + i_start;
      const double *offset = data.offsets.begin() + i_start;
      const double *y = data.Ys.begin() + i_start;

      for(uword i = 0; i < etas.n_elem; ++i, ++eta, ++weight, ++offset, ++y){
        /* Very close glm.fit implementation*/
        if(*weight <= 0.)
          continue;

        double mu = family::linkinv(*eta);
        double varmu  = family::variance(mu);
        double mu_eta_val = family::mu_eta(*eta);

        if(std::abs(mu_eta_val) < sqrt(std::numeric_limits<double>::epsilon()))
          continue;

        double z = (*eta - *offset) + (*y - mu)/mu_eta_val;
        double w = (*weight * mu_eta_val * mu_eta_val)/varmu;

        my_XtopWz += my_X.col(i) * (w * z);
        sym_mat_rank_one_update(w, my_X.col(i), my_XtopWX);
      }

      // Update shared variable
      {
        std::lock_guard<std::mutex> lk(data.m_XtopWz);
        data.XtopWz +=  my_XtopWz;
      }

      {
        std::lock_guard<std::mutex> lk(data.m_XtopWX);
        data.XtopWX += my_XtopWX;
      }
    }
  };

  static inline void compute_hessian_n_score(data_holder &data, bool first_it){
    uword n = data.X.n_cols;
    uword p = data.X.n_rows;
    data.XtopWz = arma::vec(p, arma::fill::zeros);
    data.XtopWX = arma::mat(p, p, arma::fill::zeros);

    // Compute the number of blocks to create
    uword num_blocks=(n + block_size - 1) / block_size;
    std::vector<std::future<void> > futures(num_blocks-1);
    thread_pool pool(num_blocks - 1, data.max_threads);

    std::vector<worker> workers;
    workers.reserve(num_blocks - 1);

    // declare outsite of loop to ref after loop
    uword i_start = 0;
    for(uword i = 0; i < num_blocks - 1; ++i, i_start += block_size){
      workers.emplace_back(i_start, i_start + block_size - 1, data, first_it);

      futures[i] = pool.submit(workers.back());
    }
    // compute last enteries on this thread
    worker(i_start, n - 1, data,first_it)();

    for(unsigned long i = 0; i < num_blocks - 1; ++i)
    {
      futures[i].get();   // will throw if any of the threads did
    }

    data.XtopWX = arma::symmatu(data.XtopWX);
  }

public:
  static arma::vec compute(
      arma::mat &X, arma::vec &beta0, const arma::vec &Ys,
      const arma::vec &weights, const arma::vec &offsets,
      double tol, int nthreads, int it_max, bool trace){
#if defined(USE_OPEN_BLAS)
    int openblas_nthread = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif

    /* TODO: make QR decomp of X? */

    uword p = X.n_rows;
    uword n = X.n_cols;
    data_holder data(X, Ys, weights, offsets, nthreads, p, n);

    if(p != beta0.n_elem or n != weights.n_elem or n != offsets.n_elem or n != Ys.n_elem)
      Rcpp::stop("Invalid input");

    arma::vec beta = beta0;
    int i;
    for(i = 0; i < it_max; ++i){
      arma::vec beta_old = beta;
      data.beta = &beta;

      compute_hessian_n_score(data, i == 0);

      beta = arma::solve(
        data.XtopWX, data.XtopWz, arma::solve_opts::no_approx);

      if(trace){
        Rcpp::Rcout << data.XtopWX << std::endl;
        Rcpp::Rcout << data.XtopWz << std::endl;
        Rcpp::Rcout << "it " << i << "\n"
                    << "beta_old:\t" << beta_old.t()
                    << "beta:    \t" << beta.t()
                    << "Delta norm is: " << arma::norm(beta - beta_old, 2) << std::endl;
      }

      if(arma::norm(beta - beta_old, 2) < tol) break;
    }

    if(i == it_max)
      Rcpp::stop("parallelglm did not converge");

#ifdef USE_OPEN_BLAS
    openblas_set_num_threads(openblas_nthread);
#endif

    return beta;
  }
};

namespace glm_families {

struct binomial {
  static inline double linkfun(double mu){
    return log(mu / (1 - mu));
  }

  static inline double linkinv(double eta){
    return 1 / (1 + exp(-eta));
  }

  static inline double variance(double mu){
    return mu * (1 - mu);
  }

  static inline double mu_eta(double eta){
    double exp_eta = exp(-eta);
    return (eta < -30 || eta > 30) ?
      std::numeric_limits<double>::epsilon() :
      exp_eta /((1 + exp_eta) * (1 + exp_eta));
  }

  static inline double initialize(double y, double weight){
    return linkfun((weight * y + 0.5)/(weight + 1));
  }
};

struct poisson {
  static inline double linkfun(double mu){
    return log(mu);
  }

  static inline double linkinv(double eta){
    return std::max(exp(eta), std::numeric_limits<double>::epsilon());
  }

  static inline double variance(double mu){
    return mu;
  }

  static inline double mu_eta(double eta){
    return std::max(exp(eta), std::numeric_limits<double>::epsilon());
  }

  static inline double initialize(double y, double weight){
    return linkfun(y + 0.1);
  }
};
}

// [[Rcpp::export]]
arma::vec parallelglm(
    arma::mat &X, /* Not const but will not be touched */
    const arma::vec &Ys,
    std::string family,
    arma::vec beta0, arma::vec weights, arma::vec offsets,
    double tol = 1e-8, int nthreads = 1, int it_max = 25,
    bool trace = false){
  arma::vec result;

  if(beta0.n_elem == 0)
    beta0 = arma::vec(X.n_rows, arma::fill::zeros);

  if(weights.n_elem == 0)
    weights = arma::vec(X.n_cols, arma::fill::ones);

  if(offsets.n_elem == 0)
    offsets = arma::vec(X.n_cols, arma::fill::zeros);

  if(family == "binomial"){
    result = parallelglm_class<glm_families::binomial>::compute(
      X, beta0, Ys, weights, offsets, tol, nthreads, it_max, trace);
  } else if(family == "poisson"){
    result = parallelglm_class<glm_families::poisson>::compute(
      X, beta0, Ys, weights, offsets, tol, nthreads, it_max, trace);
  } else
    Rcpp::stop("'family' not implemented");

  return result;
}
