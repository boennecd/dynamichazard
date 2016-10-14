#include "dynamichazard.h"

extern "C" {
  #include <R_ext/Constants.h> // for DOUBLE_EPS
  #include "stats/statsR.h"
  SEXP Cdqrls(SEXP x, SEXP y, SEXP tol, SEXP chk);
}

class Cdqrls_res {
public:
  Cdqrls_res(Rcpp::List res):
  coefficients(Rcpp::as<arma::vec>(res["coefficients"])),
  pivot(Rcpp::as<arma::uvec>(res["pivot"]) - 1)
  {}

  const arma::vec coefficients;
  const arma::uvec pivot; // this differs from initial values by being zero indexed
};

Cdqrls_res Cdqrls_wrapper(const arma::mat &X, const arma::vec &y,
                          const double eps, const bool check){
  Rcpp::List res = Rcpp::as<Rcpp::List>(
    Cdqrls(Rcpp::wrap(X), Rcpp::wrap(y), Rcpp::wrap(eps), Rcpp::wrap(check)));

  return(Cdqrls_res(res));
}






class dist_family {
public:
  virtual double link_func(const double&) const = 0;
  virtual double link_func_inv(const double&) const = 0;
  virtual double variance(const double&) const = 0;
  virtual double d_mu_d_eta(const double&) const = 0; // d mu / d eta
  virtual double dev_resids(const double&, const double&, const double&) const = 0;
};

class logit : public dist_family {
private:
  static constexpr double THRESH = 30.;
  static constexpr double MTHRESH = -30.;
  static constexpr double INVEPS = 1 / DOUBLE_EPS;

public:
  double link_func(const double &mu) const {
    return (log(mu / (1 - mu)));
  }

  double link_func_inv(const double &eta) const{
    double tmp = (eta < MTHRESH) ? DOUBLE_EPS :
      ((eta > THRESH) ? INVEPS : exp(eta));

    return( tmp / (1.0 + tmp));
  };

  double variance(const double &mu) const{
    return(mu * (1 - mu));
  };

  double d_mu_d_eta(const double& eta) const{
    double exp_eta = exp(eta);
    double opexp = 1 + exp_eta;

    return((eta > THRESH || eta < MTHRESH) ?  DOUBLE_EPS : exp_eta / (opexp * opexp));
  }

  double dev_resids(const double &y, const double &mu, const double &w) const{
    return (y > 0.0) ? - log(mu) : - log(1 - mu);
  };
};

class poisson : public dist_family {
public:
  double link_func(const double &mu) const{
    return log(mu);
  }

  double link_func_inv(const double &eta) const{
    return std::max(exp(eta), DOUBLE_EPS);
  }


  double variance(const double &mu) const{
    return mu;
  }

  double d_mu_d_eta(const double &eta) const{
    return std::max(exp(eta), DOUBLE_EPS);
  }

  double dev_resids(const double &y, const double &mu, const double &w) const{
    return (y > 0.0) ? w * (y * log(y / mu) - (y - mu)) : mu * w;
  }
};





using uword = arma::uword;
arma::vec IWLS(const arma::mat &X, const arma::vec &y,
               arma::vec beta, dist_family const * const family,
               const arma::vec offsets,
               const uword it_max = 1e3, const double eps = 1.0e-4){

  arma::vec weights(y.n_elem, arma::fill::ones);
  arma::vec eta = X * beta + offsets;
  arma::vec mu = eta;
  mu.for_each([&](arma::mat::value_type &val) { val = family->link_func_inv(val); });

  // TODO: delete
  // weights.for_each([&](arma::mat::value_type &val) { Rcpp::Rcout << val << "\t"; });
  // Rcpp::Rcout << std::endl;
  // eta.for_each([&](arma::mat::value_type &val) { Rcpp::Rcout << val << "\t"; });
  // Rcpp::Rcout << std::endl;
  // mu.for_each([&](arma::mat::value_type &val) { Rcpp::Rcout << val << "\t"; });
  // Rcpp::Rcout << std::endl;

  double deviance_old = 0.0;
  for(unsigned int i = 0; i < y.size(); i++)
    deviance_old += family->dev_resids(y[i], mu[i], weights[i]);

  unsigned int it = 0;
  bool converged = false;
  while(!converged && it < it_max){
    arma::vec d_mu_d_eta = eta;
    d_mu_d_eta.for_each([&](arma::vec::elem_type &val) { val = family->d_mu_d_eta(val); });

    arma::uvec good = arma::find((weights > 0) && (d_mu_d_eta != 0));
    if(good.n_elem < y.n_elem)
      Rcpp::warning("Some values 'IWLS' either had zero weights or mu = 0");

    arma::vec variance = mu.elem(good);
    variance.for_each([&](arma::vec::elem_type &val) { val = family->variance(val); });

    arma::vec z = (eta.elem(good) - offsets.elem(good)) + (y.elem(good) - mu.elem(good)) / d_mu_d_eta.elem(good);

    arma::vec w = arma::sqrt(weights.elem(good) % arma::pow(d_mu_d_eta.elem(good), 2) / variance);

    Cdqrls_res fit = Cdqrls_wrapper(
      arma::conv_to<arma::mat>::from(X.rows(good)).each_col() % w,
      z % w, eps, false);

    // TODO: deal with
    //if(any(fit$coefficients))
    // ...
    //if (nobs < fit$rank)
    // ...

    beta.elem(fit.pivot) = fit.coefficients;

    eta = X * beta + offsets;
    mu = eta;
    mu.for_each([&](arma::mat::value_type &val) { val = family->link_func_inv(val); });

    // TODO: delete
    // Rcpp::Rcout << "boh" << std::endl;
    // beta.print();
    //eta.head(10).print();

    double deviance = 0.0;
    for(unsigned int i = 0; i < y.size(); i++)
      deviance += family->dev_resids(y[i], mu[i], weights[i]);

    //TODO: Deal with
    //if (!is.finite(deviance))
    // ...
    //if (!(valideta(eta) && validmu(mu)))
    // ...

    if (std::abs(deviance - deviance_old) / (0.1 + std::abs(deviance)) < eps) {
      converged = true;
      break;
    }

    deviance_old = deviance;
    ++it;
  }

  if(!converged)
    Rcpp::stop("IWLS failed to converge");

  return beta;
}


// The export is only for internal test against glm.fit
// [[Rcpp::export]]
arma::vec IWLS_logit(const arma::mat &X, const arma::vec &y,
                     arma::vec beta,
                     const arma::vec offsets,
                     const arma::uword it_max = 1e3, const double eps = 1.0e-4){
  return(IWLS(X, y, beta, new logit(), offsets, it_max, eps));
}

// The export is only for internal test against glm.fit
// [[Rcpp::export]]
arma::vec IWLS_poisson(const arma::mat &X, const arma::vec &y,
                       arma::vec beta,
                       const arma::vec offsets,
                       const arma::uword it_max = 1e3, const double eps = 1.0e-4){
  return(IWLS(X, y, beta, new poisson(), offsets, it_max, eps));
}
