#ifndef BIGGLM_WRAPPER
#define BIGGLM_WRAPPER

#include <memory>
#include "arma_n_rcpp.h"

// Classes for R like distributions families
class dist_family {
public:
  virtual double link_func(const double&) const = 0;
  virtual double link_func_inv(const double&) const = 0;
  virtual double variance(const double&) const = 0;
  virtual double d_mu_d_eta(const double&) const = 0; // d mu / d eta
  virtual double dev_resids(const double&, const double&, const double&) const = 0;

  virtual double time_offset(const double&) const = 0;
};


// Functions and classes for fixed effects. Similar to object in bigglm
int binomialCoeff(int n, int k);

class qr_obj{
public:
  qr_obj(unsigned int p):
  D(new arma::vec(p, arma::fill::zeros)), rbar(new arma::vec((p == 1)? 0 : binomialCoeff(p, 2), arma::fill::zeros)),
    thetab(new arma::vec(p, arma::fill::zeros)), ss(0.), checked(false),
    tol(new arma::vec(p, arma::fill::zeros))
  {}
  qr_obj() = default;

  std::shared_ptr<arma::vec> D;
  std::shared_ptr<arma::vec> rbar;
  std::shared_ptr<arma::vec> thetab;
  double ss;
  bool checked;
  std::shared_ptr<arma::vec> tol;
};


template<class T>
class bigglm_updateQR{
  // match logic in update.bigqr
  arma::vec linkinv(const arma::vec &eta);

  arma::vec d_mu_d_eta(const arma::vec &eta);

  arma::vec variance(const arma::vec &mu);

protected:
  T t;

public:
  bigglm_updateQR<T>(): t() {}

  void update(qr_obj &qr, // Previous/starting value. Will be overwritten
              const arma::mat &X, const arma::vec &eta,
              const arma::vec &offset, arma::vec &y, // y will not be altered
              const arma::vec &w);
};

arma::vec bigglm_regcf(qr_obj &qr);

//logit_fam
class logit_fam : public dist_family {
private:
  static constexpr double THRESH = 30.;
  static constexpr double MTHRESH = -30.;
  static constexpr double INVEPS = 1 / DOUBLE_EPS;

public:
  double link_func(const double &mu) const{
    return (log(mu / (1 - mu)));
  };

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
  };

  double dev_resids(const double &y, const double &mu, const double &w) const{
    return (y > 0.0) ? - log(mu) : - log(1 - mu);
  };

  double time_offset(const double &delta_t) const{
    return 0.;
  }
};

//poisson_fam
class poisson_fam : public dist_family
{
public:
  double link_func(const double &mu) const{
    return log(mu);
  };

  double link_func_inv(const double &eta) const{
    return std::max(exp(eta), DOUBLE_EPS);
  };

  double variance(const double &mu) const{
    return mu;
  };

  double d_mu_d_eta(const double &eta) const{
    return std::max(exp(eta), DOUBLE_EPS);
  };

  double dev_resids(const double &y, const double &mu, const double &w) const{
    return (y > 0.0) ? w * (y * log(y / mu) - (y - mu)) : mu * w;
  };

  double time_offset(const double &delta_t) const{
    return log(delta_t);
  };
};


// Define the concrete templates
template class bigglm_updateQR<logit_fam>;
template class bigglm_updateQR<poisson_fam>;

#endif
