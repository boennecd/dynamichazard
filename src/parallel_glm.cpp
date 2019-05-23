/* Different methods to fit a GLM model in parallel. Both are much like the
 * `bam` method in the `mgcv` package.
 *
 * TODO: there is some redundancy in the parallel_glm classes. Make a base
 *       class which the other two derives from?                             */

#include "arma_n_rcpp.h"
#include "thread_pool.h"
#include "parallel_qr.h"
#include "arma_BLAS_LAPACK.h"
#include "family.h"

#define COMPUTE_ARGS                                           \
  arma::mat &X,                                                \
  arma::vec &beta0,                                            \
  arma::vec &Ys,                                               \
  arma::vec &weights,                                          \
  arma::vec &offsets,                                          \
  std::unique_ptr<glm_base> family,                            \
  double tol,                                                  \
  int nthreads,                                                \
  int it_max,                                                  \
  bool trace

#define MIN(a,b) (((a)<(b))?(a):(b))

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
  const std::unique_ptr<glm_base> family;
  const arma::uword block_size;

  data_holder_base(
    arma::mat &X, arma::vec &Ys, arma::vec &weights, arma::vec &offsets,
    const arma::uword max_threads, const arma::uword p, const arma::uword n,
    std::unique_ptr<glm_base> family, arma::uword block_size = 10000):
    X(X), Ys(Ys), weights(weights), offsets(offsets),
    max_threads(max_threads), p(p), n(n), family(std::move(family)),
    block_size(block_size)
  {}
};


struct parallelglm_out {
  arma::vec result;
  unsigned int iter;
};

/* Similar function glm.fit in r-source/src/library/stats/R/glm.R
 * It only returns the coefficient vector. Computes X^\top W X and may be less
 * stable */
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

  /* worker class for multithreading */
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
          *eta = data.family->glm_initialize(*y, *weight);

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

        double mu = data.family->glm_linkinv(*eta);
        double varmu  = data.family->glm_variance(mu);
        double mu_eta_val = data.family->glm_mu_eta(*eta);

        if(std::abs(mu_eta_val) < std::numeric_limits<double>::epsilon())
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
  static parallelglm_out compute(COMPUTE_ARGS){
    uword p = X.n_rows;
    uword n = X.n_cols;
    data_holder data(
        X, Ys, weights, offsets, nthreads, p, n, std::move(family));

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

    parallelglm_out out;
    out.result = beta;
    out.iter = (unsigned int)MIN(i + 1L, it_max);

    return out;
  }
};

/* Class to fit glm using QR updated in chunks */
class parallelglm_class_QR {
  using uword = arma::uword;

  class glm_qr_data_generator : public qr_data_generator {
    static constexpr double zero_eps = 1e-100;

    const uword i_start, i_end;
    data_holder_base &data;
    const bool first_it;

  public:
    glm_qr_data_generator
    (uword i_start, uword i_end, data_holder_base &data, bool first_it):
    i_start(i_start), i_end(i_end), data(data), first_it(first_it) {}

    qr_work_chunk get_chunk() const override {
      /* assign objects for later use */
      arma::span my_span(i_start, i_end);
      uword n = i_end - i_start + 1;

      arma::vec y     (data.Ys.begin()      + i_start           , n, false);
      arma::vec weight(data.weights.begin() + i_start           , n, false);
      arma::vec offset(data.offsets.begin() + i_start           , n, false);
      arma::mat X     (data.X.begin() + i_start * data.p, data.p, n);

      arma::vec eta;
      if(first_it){
        eta = arma::vec(i_end - i_start + 1);
        double *eta_i = eta.begin();
        const double *y_i = y.begin();
        const double *wt = weight.begin();
        for(uword i = 0; i < n; ++i, ++eta_i, ++wt, ++y_i)
          *eta_i = data.family->glm_initialize(*y_i, *wt);

      } else
        eta = (data.beta->t() * X).t() + offset;

      /* compute values for QR computation */
      arma::vec mu = eta, mu_eta_val = eta;
      double *m = mu.begin(), *mev = mu_eta_val.begin();
      for(uword i = 0; i < mu.n_elem; ++i, ++m, ++mev){
        *m = data.family->glm_linkinv(*m);
        *mev = data.family->glm_mu_eta(*mev);
      }

      arma::uvec good = arma::find(
        (weight > 0) %
          ((-zero_eps < mu_eta_val) + (mu_eta_val < zero_eps) != 2));

      const double *mu_i = mu.begin();
      const double *wt_i = weight.begin();
      const double  *y_i = y.begin();

      double dev = 0;
      for(uword i = 0; i < n; ++i, ++mu_i, ++wt_i, ++y_i)
        dev += data.family->glm_dev_resids(*y_i, *mu_i, *wt_i);

      mu = mu(good);
      eta = eta(good);
      mu_eta_val = mu_eta_val(good);
      arma::vec var = mu;
      for(auto v = var.begin(); v != var.end(); ++v)
        *v = data.family->glm_variance(*v);

      /* compute X and working responses and return */
      arma::vec z = (eta - offset(good)) + (y(good) - mu) / mu_eta_val;
      arma::vec w = arma::sqrt(
        (weight(good) % arma::square(mu_eta_val)) / var);

      X = X.cols(good);
      X = X.t();
      X.each_col() %= w;
      z %= w;

      arma::mat dev_mat; dev_mat = dev; /* make 1 x 1 matrix */

      return { std::move(X), std::move(z), dev_mat};
    }
  };

public:
  static parallelglm_out compute(COMPUTE_ARGS){
    uword p = X.n_rows;
    uword n = X.n_cols;
    data_holder_base data(
        X, Ys, weights, offsets, nthreads, p, n, std::move(family));

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

      R_F R_f_out = get_R_f(data, i == 0);
      /* TODO: can maybe done smarter using that R is triangular before
       *       permutation */
      arma::mat R = R_f_out.R_rev_piv();
      beta = arma::solve(R.t(), R.t() * R_f_out.F.col(0),
                         arma::solve_opts::no_approx);
      beta = arma::solve(R    , beta,
                         arma::solve_opts::no_approx);

      if(trace){
        Rcpp::Rcout << "it " << i << "\n"
                    << "beta_old:\t" << beta_old.t()
                    << "beta:    \t" << beta.t()
                    << "Delta norm is: "<< std::endl
                    << arma::norm(beta - beta_old, 2) << std::endl
                    << "deviance is " << dev << std::endl;
      }

      double devold = dev;
      dev = R_f_out.dev(0, 0);

      if(std::abs(dev - devold) / (1e-2 + std::abs(dev)) < tol)
        break;
    }

    if(i == it_max)
      Rcpp::stop("parallelglm did not converge");

    parallelglm_out out;
    out.result = beta;
    out.iter = (unsigned int)MIN(i + 1L, it_max);

    return out;
  }

  static R_F get_R_f(data_holder_base &data, bool first_it){
    uword n = data.X.n_cols;

    // Compute the number of blocks to create
    uword num_blocks = (n + data.block_size - 1) / data.block_size;
    std::vector<std::unique_ptr<qr_data_generator>> generators;
    generators.reserve(num_blocks);

    // setup generators
    uword i_start = 0;
    for(uword i = 0; i < num_blocks; ++i, i_start += data.block_size)
      generators.emplace_back(
        new glm_qr_data_generator(
              i_start, std::min(n - 1, i_start + data.block_size - 1),
              data, first_it));

    return qr_parallel(std::move(generators), data.max_threads).compute();
  }
};


template<class TMethod>
parallelglm_out parallelglm_fit(
    arma::mat &X, arma::vec &Ys, std::string family, arma::vec beta0,
    arma::vec &weights, arma::vec &offsets, double tol, int nthreads,
    int it_max, bool trace){
  parallelglm_out result;
  if(family == BINOMIAL)
    result = TMethod::compute(
      X, beta0, Ys, weights, offsets,
      std::unique_ptr<glm_base>(new logistic()), tol, nthreads, it_max, trace);
  else if(family == CLOGLOG)
    result = TMethod::compute(
      X, beta0, Ys, weights, offsets,
      std::unique_ptr<glm_base>(new cloglog()), tol, nthreads, it_max, trace);
  else if(family == POISSON)
    result = TMethod::compute(
      X, beta0, Ys, weights, offsets,
      std::unique_ptr<glm_base>(new exponential()), tol, nthreads, it_max,
      trace);
  else
    Rcpp::stop("'family' not implemented");

  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector parallelglm(
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

  parallelglm_out result;
  if(method == "quick")
    result = parallelglm_fit<parallelglm_class_quick>(
      X, Ys, family, beta0, weights, offsets, tol, nthreads, it_max, trace);
  else if(method == "QR")
    result = parallelglm_fit<parallelglm_class_QR>(
      X, Ys, family, beta0, weights, offsets, tol, nthreads, it_max, trace);
  else
    Rcpp::stop("'method' not implemented");

  Rcpp::NumericVector out =
    Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(result.result));
  out.attr("iter") = result.iter;

  return out;
}
