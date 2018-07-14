#include "arma_n_rcpp.h"
#include "parallel_qr.h"
#include "arma_BLAS_LAPACK.h"
#include "family.h"

/* see https://stackoverflow.com/a/18776112/5861244 */
inline void
  copy
  (double *out, const double *org, const unsigned int org_n_elem,
   unsigned int n_times){
    unsigned int num_copied = org_n_elem, num_total = org_n_elem * n_times;
    memcpy(out, org, num_copied * sizeof(double));

    while(num_copied * 2 <= num_total) {
      memcpy(out + num_copied, out, num_copied * sizeof(double));
      num_copied *= 2;
    }

    if(num_copied < num_total)
      memcpy(out + num_copied, out,
             (num_total - num_copied) * sizeof(double));
  }

inline arma::mat arma_copy(const arma::mat &org, const arma::uword n_times){
  arma::mat out(org.n_rows, org.n_cols * n_times);

  copy(&out[0], &org[0], org.n_elem, n_times);

  return out;
}

inline arma::vec arma_copy(const arma::vec &org, const arma::uword n_times){
  arma::vec out(org.n_elem * n_times);

  copy(&out[0], &org[0], org.n_elem, n_times);

  return out;
}

inline arma::vec set_eta
  (const arma::mat &X, const arma::vec &beta, const arma::uword n_times){
  arma::vec tmp = X.t() * beta;
  return arma_copy(tmp, n_times);
}

inline arma::vec set_offsets
  (const arma::vec &dts, const glm_base *family, const arma::uword n_times){
  arma::vec offsets;

  if(family->name() == POISSON){
    offsets = arma::log(dts);
    offsets = arma_copy(offsets, n_times);

  } else if(family->name() == BINOMIAL){
    offsets = arma::vec(dts.n_elem * n_times, arma::fill::zeros);

  } else
    Rcpp::stop("family not implemented");

  return offsets;
}

struct data_holder {
  const arma::uword n_obs;
  const arma::vec eta_wo_offset;
  const arma::vec at_risk_offsets;
  const arma::mat X;
  arma::vec Y;   /* only non-const to use non-copy constructor later */
  arma::vec dts; /* only non-const to use non-copy constructor later */
  const arma::mat &cloud;
  const arma::vec &cl_weights;
  const arma::mat &ran_vars;
  const arma::vec &beta;
  std::unique_ptr<glm_base> family;
  const arma::uword n_times;

  data_holder
    (const arma::mat &X, const arma::vec &Y, const arma::vec &dts,
     const arma::mat &cloud, const arma::vec &cl_weights,
     const arma::mat &ran_vars, const arma::vec &beta,
     std::unique_ptr<glm_base> family, const arma::uword n_times):
    n_obs(X.n_cols), eta_wo_offset(set_eta(X, beta, n_times)),
    at_risk_offsets(set_offsets(dts, family.get(), n_times)),
    X(arma_copy(X, n_times)), Y(arma_copy(Y, n_times)),
    dts(arma_copy(dts, n_times)), cloud(cloud), cl_weights(cl_weights),
    ran_vars(ran_vars), beta(beta), family(std::move(family)), n_times(n_times)
    { }
};

class pf_fixed_generator : public qr_data_generator {
  using uword = arma::uword;
  static constexpr double zero_eps = 1e-100;

  data_holder &data; /* non-const to use non-copy constructor later for arma objects */
  const uword istart, iend;

public:
  pf_fixed_generator
  (data_holder &data, const uword istart, const uword iend):
  data(data), istart(istart), iend(iend) { }

  qr_work_chunk get_chunk() const override {
    /* assign objects for later use */
    arma::mat X;
    arma::vec Y, dts, eta, offset;
    arma::uword n_particles = iend - istart + 1L;

    if(n_particles < data.n_times){
      arma::uword n_keep = n_particles * data.n_obs;

      X = arma::mat(data.X.begin(), data.X.n_rows, n_keep);
      eta = arma::vec(data.eta_wo_offset.begin(), n_keep);
      offset = arma::vec(data.at_risk_offsets.begin(), n_keep);
      Y   = arma::vec(data.Y.begin()  , n_keep, false);
      dts = arma::vec(data.dts.begin(), n_keep, false);

    } else {
      X = data.X;                    // copy
      eta = data.eta_wo_offset;      // copy
      offset = data.at_risk_offsets; // copy
      Y   = arma::vec(data.Y.begin()  , data.Y.n_elem  , false);
      dts = arma::vec(data.dts.begin(), data.dts.n_elem, false);

    }
    uword n = X.n_cols;
    arma::vec weight(n);

    for(arma::uword i = 0; i < n_particles; ++i){
      weight.subvec(i * data.n_obs, (i + 1L) * data.n_obs - 1L).fill(
          data.cl_weights[istart + i]);
      offset.subvec(i * data.n_obs, (i + 1L) * data.n_obs - 1L) +=
        data.ran_vars.t() * data.cloud.col(istart + i);

    }

    /* compute values for QR computation */
    eta += offset;
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
    const double  *y_i = Y.begin();

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
    arma::vec z = (eta - offset(good)) + (data.Y(good) - mu) / mu_eta_val;
    arma::vec w = arma::sqrt(
      (weight(good) % arma::square(mu_eta_val)) / var);

    X = X.cols(good);
    arma::inplace_trans(X);
    X.each_col() %= w;
    z %= w;

    arma::mat dev_mat(1, 1);
    dev_mat(0,0) = dev;

    return { std::move(X), std::move(z), dev_mat };
  }
};


// [[Rcpp::export]]
Rcpp::NumericMatrix test_copy_mat(const arma::mat &X, const int n_times){
  return(Rcpp::wrap(arma_copy(X, n_times)));
}

// [[Rcpp::export]]
Rcpp::NumericVector test_copy_vec(const arma::vec &x, const int n_times){
  return(Rcpp::wrap(arma_copy(x, n_times)));
}

// [[Rcpp::export]]
Rcpp::List pf_fixed_effect_iteration(
    const arma::mat &X, const arma::vec &Y, const arma::vec &dts,
    const arma::mat &cloud, const arma::vec &cl_weights,
    const arma::mat &ran_vars, const arma::vec &beta,
    std::string family, int max_threads, const bool debug,
    const unsigned int max_bytes = 1000000){
  arma::uword n_particles = cloud.n_cols;
  arma::uword max_blocks = std::min(
    std::max((long int)(max_bytes / (X.n_rows * X.n_cols * 8L)), 1L),
    (long int)n_particles);
  arma::uword even_blocks = n_particles / max_threads + 1L;
  arma::uword block_use = std::min(max_blocks, even_blocks);

  if(debug)
    Rcpp::Rcout << "Making chunks with " << std::setw(12)
                << block_use << " particles out of "
                << std::setw(12) << n_particles << " particles in "
                << "`pf_fixed_effect_iteration`" << std::endl;

  std::unique_ptr<data_holder> dat;
  if(family == BINOMIAL){
    dat.reset(new data_holder(
        X, Y, dts, cloud, cl_weights, ran_vars, beta,
        std::unique_ptr<glm_base>(new logistic()), block_use));

  } else if(family == POISSON){
    dat.reset(new data_holder(
        X, Y, dts, cloud, cl_weights, ran_vars, beta,
        std::unique_ptr<glm_base>(new exponential()), block_use));

  } else
    Rcpp::stop("Family not implemented");

  /* setup generators */
  std::vector<std::unique_ptr<qr_data_generator>> generators;
  for(arma::uword i_start = 0L; i_start < n_particles; i_start += block_use)
    generators.emplace_back(new pf_fixed_generator(
        *dat.get(), i_start, std::min(
            n_particles - 1L, i_start + block_use - 1L)));

  // compute and return
  R_F res = qr_parallel(std::move(generators), max_threads).compute();

  return Rcpp::List::create(
    Rcpp::Named("R") =  Rcpp::wrap(res.R),
    Rcpp::Named("pivot") =  Rcpp::wrap(res.pivot),
    Rcpp::Named("f") =  Rcpp::wrap(res.F),
    Rcpp::Named("dev") =  Rcpp::wrap(res.dev));
}
