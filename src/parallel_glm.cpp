/* Different methods to fit a GLM model in parallel. Both are much like the
 * `bam` method in the `mgcv` package.
 *
 * TODO: there is some redundancy in the parallel_glm classes. Make a base
 *       class which the other two derives from.                             */

#include "arma_n_rcpp.h"
#include "thread_pool.h"
#include "arma_BLAS_LAPACK.h"

#define COMPUTE_ARGS                                           \
  arma::mat &X,                                                \
  arma::vec &beta0,                                            \
  arma::vec &Ys,                                               \
  arma::vec &weights,                                          \
  arma::vec &offsets,                                          \
  double tol,                                                  \
  int nthreads,                                                \
  int it_max,                                                  \
  bool trace

#define BINOMIAL "binomial"
#define POISSON "poisson"

/* base data holder class */
/* data holder class */
class data_holder_base {
public:
  arma::vec *beta;

  /* These are not const but should not be changed... */
  arma::mat &X;
  arma::vec &Ys;
  arma::vec &weights;
  arma::vec &offsets;

  const arma::uword max_threads, p, n;
  const arma::uword block_size;

  data_holder_base(
    arma::mat &X, arma::vec &Ys, arma::vec &weights, arma::vec &offsets,
    const arma::uword max_threads, const arma::uword p, const arma::uword n,
    arma::uword block_size = 10000):
    X(X), Ys(Ys), weights(weights), offsets(offsets),
    max_threads(max_threads), p(p), n(n), block_size(block_size)
  {}
};

/* Similar function glm.fit in r-source/src/library/stats/R/glm.R
 * It only returns the coefficient vector. Computes X^\top W X and may be less
 * stable */

template<typename family>
class parallelglm_class_quick {
  using uword = arma::uword;

  class data_holder : public data_holder_base {
  public:
    arma::mat XtopWX;
    arma::vec XtopWz;
    std::mutex m_XtopWX;
    std::mutex m_XtopWz;

    using data_holder_base::data_holder_base;
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

  static void compute_hessian_n_score(data_holder &data, bool first_it){
    uword n = data.X.n_cols;
    uword p = data.X.n_rows;
    data.XtopWz = arma::vec(p, arma::fill::zeros);
    data.XtopWX = arma::mat(p, p, arma::fill::zeros);

    // Compute the number of blocks to create
    uword num_blocks=(n + data.block_size - 1) / data.block_size;
    std::vector<std::future<void> > futures(num_blocks-1);
    thread_pool pool(std::min(num_blocks - 1, data.max_threads));

    std::vector<worker> workers;
    workers.reserve(num_blocks - 1);

    // declare outsite of loop to ref after loop
    uword i_start = 0;
    for(uword i = 0; i < num_blocks - 1; ++i, i_start += data.block_size){
      workers.emplace_back(i_start, i_start + data.block_size - 1, data, first_it);

      futures[i] = pool.submit(workers.back());
    }
    // compute last enteries on this thread
    worker(i_start, n - 1, data, first_it)();

    for(unsigned long i = 0; i < num_blocks - 1; ++i)
    {
      futures[i].get();   // will throw if any of the threads did
    }

    data.XtopWX = arma::symmatu(data.XtopWX);
  }

public:
  static arma::vec compute(COMPUTE_ARGS){
    uword p = X.n_rows;
    uword n = X.n_cols;
    data_holder data(X, Ys, weights, offsets, nthreads, p, n);

    if(p != beta0.n_elem or
         n != weights.n_elem or
         n != offsets.n_elem or n != Ys.n_elem)
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
                    << "Delta norm is: "
                    << arma::norm(beta - beta_old, 2) << std::endl;
      }

      if(arma::norm(beta - beta_old, 2) < tol) break;
    }

    if(i == it_max)
      Rcpp::stop("parallelglm did not converge");

    return beta;
  }
};

/* Class to fit glm using QR updated in chunks */

/* Let x = QR. Then this is a data holder for the R matrix and f = Q^\top y.
 * dev is deviance computed for the chunk computed with the current
 * coefficient vector                                                        */
struct R_f {
  const arma::mat R;
  const arma::uvec pivot;
  const arma::vec f;
  const double dev;

  arma::mat R_rev_piv() const {
    arma::uvec piv = pivot;
    piv(piv) = arma::regspace<arma::uvec>(0, 1, piv.n_elem - 1);
    return R.cols(piv);
  }
};

/* function to combine results */
R_f R_f_combine(std::vector<R_f> &R_fs){
  R_f *R_fs_i = &R_fs.back();
  arma::mat R_stack = R_fs_i->R_rev_piv();
  arma::vec f_stack = std::move(R_fs_i->f);
  double  dev = R_fs_i->dev;
  R_fs.pop_back();

  while(!R_fs.empty()){
    R_fs_i = &R_fs.back();
    R_stack = arma::join_cols(R_stack, R_fs_i->R_rev_piv());
    f_stack = arma::join_cols(f_stack, std::move(R_fs_i->f));
    dev += R_fs_i->dev;

    R_fs.pop_back();
  }

  /* make new QR decomp and compute new f*/
  QR_factorization qr(R_stack);
  arma::vec f = qr.qy(f_stack, true).subvec(0, R_stack.n_cols - 1);

  return R_f { qr.R(), qr.pivot(), std::move(f), dev };
}

template<typename family>
class parallelglm_class_QR {
  using uword = arma::uword;

  /* worker class for multithreading*/
  class worker {
    static constexpr double zero_eps = 2.220446e-16;

    const uword i_start, i_end;
    data_holder_base &data;
    const bool first_it;

  public:
    worker(uword i_start, uword i_end, data_holder_base &data, bool first_it):
    i_start(i_start), i_end(i_end), data(data), first_it(first_it) {}

    R_f operator()(){

      /* assign objects for later use */
      arma::span my_span(i_start, i_end);
      uword n = i_end - i_start + 1;

      arma::vec y     (data.Ys.begin()      + i_start, n, false);
      arma::vec weight(data.weights.begin() + i_start, n, false);
      arma::vec offset(data.offsets.begin() + i_start, n, false);
      arma::mat X     (data.X.begin() + i_start * data.p, data.p, n);

      arma::vec eta;
      if(first_it){
        eta = arma::vec(i_end - i_start + 1);
        double *eta_i = eta.begin();
        const double *y_i = y.begin();
        const double *wt = weight.begin();
        for(uword i = 0; i < n; ++i, ++eta_i, ++wt, ++y_i)
          *eta_i = family::initialize(*y_i, *wt);

      } else
        eta = (data.beta->t() * X).t() + offset;

      /* compute values for QR computation */
      arma::vec mu = eta, mu_eta_val = eta;
      mu.transform(family::linkinv);
      mu_eta_val.transform(family::mu_eta);

      arma::uvec good = arma::find(
        (weight > 0) %
          ((-zero_eps < mu_eta_val) + (mu_eta_val < zero_eps) != 2));

      /* compute deviance */
      const double *mu_i = mu.begin();
      const double *wt_i = weight.begin();
      const double  *y_i = y.begin();

      double dev = 0;
      for(uword i = 0; i < n; ++i, ++mu_i, ++wt_i, ++y_i)
        dev += family::dev_resids(*y_i, *mu_i, *wt_i);

      /* make QR decomposition */
      mu = mu(good);
      eta = eta(good);
      mu_eta_val = mu_eta_val(good);
      arma::vec var = mu;
      var.transform(family::variance);

      arma::vec z = (eta - offset(good)) + (y(good) - mu) / mu_eta_val;
      arma::vec w = arma::sqrt(
        (weight(good) % arma::square(mu_eta_val)) / var);

      /* find QR */
      X = X.cols(good);
      X = X.t();
      X.each_col() %= w;
      QR_factorization qr(X);

      z %= w;
      arma::vec f = qr.qy(z, true).subvec(0, data.p - 1);

      return R_f { qr.R(), qr.pivot(), std::move(f), dev };
    }
  };

public:
  static arma::vec compute(COMPUTE_ARGS){
    uword p = X.n_rows;
    uword n = X.n_cols;
    data_holder_base data(X, Ys, weights, offsets, nthreads, p, n);

    if(p != beta0.n_elem or
         n != weights.n_elem or
         n != offsets.n_elem or n != Ys.n_elem)
      Rcpp::stop("Invalid input");

    arma::vec beta = beta0;
    int i;
    double dev = 0;
    for(i = 0; i < it_max; ++i){
      arma::vec beta_old = beta;
      data.beta = &beta;

      R_f R_f_out = get_R_f(data, i == 0);
      arma::mat R = R_f_out.R_rev_piv();
      beta = arma::solve(R.t() * R, R.t() * R_f_out.f);

      if(trace){
        Rcpp::Rcout << "it " << i << "\n"
                    << "beta_old:\t" << beta_old.t()
                    << "beta:    \t" << beta.t()
                    << "Delta norm is: "<< std::endl
                    << arma::norm(beta - beta_old, 2) << std::endl
                    << "deviance is " << dev << std::endl;
      }

      double devold = dev;
      dev = R_f_out.dev;

      if(std::abs(dev - devold) / (1e-2 + std::abs(dev)) < tol)
        break;
    }

    if(i == it_max)
      Rcpp::stop("parallelglm did not converge");

    return beta;
  }

  static R_f get_R_f(data_holder_base &data, bool first_it){
    uword n = data.X.n_cols;

    // Compute the number of blocks to create
    uword num_blocks=(n + data.block_size - 1) / data.block_size;
    std::vector<std::future<R_f> > futures(num_blocks-1);
    thread_pool pool(std::min(num_blocks - 1, data.max_threads));

    std::vector<worker> workers;
    workers.reserve(num_blocks - 1);

    // declare outsite of loop to ref after loop
    uword i_start = 0;
    for(uword i = 0; i < num_blocks - 1; ++i, i_start += data.block_size){
      workers.emplace_back(i_start, i_start + data.block_size - 1, data, first_it);

      futures[i] = pool.submit(workers.back());
    }

    // compute last enteries on this thread
    std::vector<R_f> R_fs;
    R_fs.push_back(worker(i_start, n - 1, data, first_it)());

    for(unsigned long i = 0; i < num_blocks - 1; ++i)
    {
      R_fs.push_back(futures[i].get());
    }

    // Find final R, f and deviance
    return R_f_combine(R_fs);
  }
};



/* glm families */
namespace glm_families {
  struct binomial {
    static inline double dev_resids(double y, double mu, double wt){
      return - 2 * wt * (y * log(mu) + (1 - y) * log(1 - mu));
    }

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
    static inline double dev_resids(double y, double mu, double wt){
      if(y > 0)
        return 2 * (wt * (y * log(y/mu) - (y - mu)));

      return 2 * mu * wt;
    }

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

template<template<typename> class TMethod>
arma::vec parallelglm_fit(
    arma::mat &X, arma::vec &Ys, std::string family, arma::vec beta0,
    arma::vec &weights, arma::vec &offsets, double tol, int nthreads,
    int it_max, bool trace){
  arma::vec result;
  if(family == BINOMIAL)
    result = TMethod<glm_families::binomial>::compute(
      X, beta0, Ys, weights, offsets, tol, nthreads, it_max, trace);
  else if(family == POISSON)
    result = TMethod<glm_families::poisson>::compute(
      X, beta0, Ys, weights, offsets, tol, nthreads, it_max, trace);
  else
    Rcpp::stop("'family' not implemented");

  return result;
}

// [[Rcpp::export]]
arma::vec parallelglm(
    arma::mat &X, arma::vec &Ys, std::string family, arma::vec beta0,
    arma::vec &weights, arma::vec &offsets, double tol = 1e-8,
    int nthreads = 1, int it_max = 25, bool trace = false,
    std::string method = "Quick"){
  if(beta0.n_elem == 0)
    beta0 = arma::vec(X.n_rows, arma::fill::zeros);

  if(weights.n_elem == 0)
    weights = arma::vec(X.n_cols, arma::fill::ones);

  if(offsets.n_elem == 0)
    offsets = arma::vec(X.n_cols, arma::fill::zeros);

  arma::vec result;
  if(method == "quick")
    result = parallelglm_fit<parallelglm_class_quick>(
      X, Ys, family, beta0, weights, offsets, tol, nthreads, it_max, trace);
  else if(method == "QR")
    result = parallelglm_fit<parallelglm_class_QR>(
      X, Ys, family, beta0, weights, offsets, tol, nthreads, it_max, trace);
  else
    Rcpp::stop("'method' not implemented");

  return result;
}

// Exported to test the intermediate computations
// [[Rcpp::export]]
Rcpp::List parallelglm_QR_test(
    arma::mat &X, arma::vec &Ys, std::string family, arma::vec beta0,
    arma::vec &weights, arma::vec &offsets, double tol = 1e-8,
    int nthreads = 1, int it_max = 25, bool trace = false,
    int block_size = 100){
  arma::uword p = X.n_rows;
  arma::uword n = X.n_cols;
  data_holder_base data(X, Ys, weights, offsets, nthreads, p, n, block_size);
  data.beta = &beta0;

  std::unique_ptr<R_f> res;
  if(family == BINOMIAL){
    res.reset(new R_f(parallelglm_class_QR<glm_families::binomial>::get_R_f(
      data, false /* note the false */)));
  } else if(family == POISSON){
    res.reset(new R_f(parallelglm_class_QR<glm_families::poisson>::get_R_f(
      data, false /* note the false */)));
  }

  return Rcpp::List::create(
    Rcpp::Named("R") =  Rcpp::wrap(res->R),
    Rcpp::Named("pivot") =  Rcpp::wrap(res->pivot),
    Rcpp::Named("f") =  Rcpp::wrap(res->f),
    Rcpp::Named("dev") =  Rcpp::wrap(res->dev));
}
