/*
 Similar function glm.fit in r-source/src/library/stats/R/glm.R
 It only returns the coefficient vector
*/

#include "arma_n_rcpp.h"
#include "thread_pool.h"
#include "arma_utils.h"

template<typename family>
class parallelglm_class {
  using uword = arma::uword;
  static const uword block_size = 1000;

  /* data holder class */
  struct data_holder {
    const arma::mat &X;
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
      const arma::mat &X, const arma::vec &Ys, const arma::vec &weights,
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

  public:
    worker(uword i_start, uword i_end, data_holder &data):
    i_start(i_start), i_end(i_end), data(data) {}

    void operator()(){

      arma::span my_span(i_start, i_end);
      arma::mat my_X = data.X(arma::span::all, my_span);
      arma::vec etas = (data.beta->t() * my_X).t() + data.offsets(my_span);

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

        if(mu_eta_val == 0.)
          continue;

        double z = (*eta - *offset) + (*y - mu)/mu_eta_val;
        double w = (*weight * pow(mu_eta_val, 2))/varmu;

        /* TODO: delete
        Rcpp::Rcout << "eta "<< *eta << "\t offset " << *offset << "\t g "<< mu << "\t gprime "<< mu_eta_val << "\t var "<< varmu << "\t z "<<
          z << "\t w " << w << "\t" << std::endl;
        */

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

  static inline void compute_hessian_n_score(data_holder &data){
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
      workers.emplace_back(i_start, i_start + block_size - 1, data);

      futures[i] = pool.submit(workers.back());
    }
    // compute last enteries on this thread
    worker(i_start, n - 1, data)();

    for(unsigned long i = 0; i < num_blocks - 1; ++i)
    {
      futures[i].get();   // will throw if any of the threads did
    }

    data.XtopWX = arma::symmatu(data.XtopWX);
  }

public:
  static arma::vec compute(
      const arma::mat &X, arma::vec &beta0, const arma::vec &Ys,
      const arma::vec &weights, const arma::vec &offsets,
      double tol, int nthreads, int it_max){
#if defined(USE_OPEN_BLAS)
    int openblas_nthread = openblas_get_num_threads();
    openblas_set_num_threads(1);
#endif

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

      compute_hessian_n_score(data);

      beta = arma::solve(
        data.XtopWX, data.XtopWz, arma::solve_opts::no_approx);

      if(sqrt(arma::norm(beta - beta_old, 2)) < tol) break;
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
  struct binomial{
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
  };
}

// [[Rcpp::export]]
arma::vec parallelglm(
    const arma::mat &X, const arma::vec &Ys,
    std::string family,
    arma::vec beta0, arma::vec weights, arma::vec offsets,
    double tol = 1e-8, int nthreads = 1, int it_max = 25){
  arma::vec result;

  if(beta0.n_elem == 0)
    beta0 = arma::vec(X.n_rows, arma::fill::zeros);

  if(weights.n_elem == 0)
    weights = arma::vec(X.n_cols, arma::fill::ones);

  if(offsets.n_elem == 0)
    offsets = arma::vec(X.n_cols, arma::fill::zeros);

  if(family == "binomial"){
    result = parallelglm_class<glm_families::binomial>::compute(
      X, beta0, Ys, weights, offsets, tol, nthreads, it_max);
  } else
    Rcpp::stop("'family' not implemented");

  return result;
}
