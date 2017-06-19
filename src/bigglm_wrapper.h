#ifndef BIGGLM_WRAPPER
#define BIGGLM_WRAPPER

#include <memory>
#include "arma_n_rcpp.h"

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


template<typename T>
class bigglm_updateQR{
  // match logic in update.bigqr
  static arma::vec linkinv(const arma::vec &eta);

  static arma::vec d_mu_d_eta(const arma::vec &eta);

  static arma::vec variance(const arma::vec &mu);

public:
  static void update(
      qr_obj &qr, // Previous/starting value. Will be overwritten
      const arma::mat &X, const arma::vec &eta,
      const arma::vec &offset, arma::vec &y, // y will not be altered
      const arma::vec &w);

  typedef T family;
};

arma::vec bigglm_regcf(qr_obj &qr);

//logit_fam
class logit_fam {
private:
  static constexpr double THRESH = 30.;
  static constexpr double MTHRESH = -30.;
  static constexpr double INVEPS = 1 / DOUBLE_EPS;

public:
  static inline double link_func(const double &mu){
    return (log(mu / (1 - mu)));
  };

  static inline double link_func_inv(const double &eta){
    double tmp = (eta < MTHRESH) ? DOUBLE_EPS :
    ((eta > THRESH) ? INVEPS : exp(eta));

    return( tmp / (1.0 + tmp));
  };

  static inline double variance(const double &mu){
    return(mu * (1 - mu));
  };

  static inline double d_mu_d_eta(const double& eta){
    double exp_eta = exp(eta);
    double opexp = 1 + exp_eta;

    return((eta > THRESH || eta < MTHRESH) ?  DOUBLE_EPS : exp_eta / (opexp * opexp));
  };

  static inline double dev_resids(const double &y, const double &mu, const double &w){
    return (y > 0.0) ? - log(mu) : - log(1 - mu);
  };

  static inline double time_offset(const double &delta_t){
    return 0.;
  }
};

//poisson_fam
class poisson_fam
{
public:
  static inline double link_func(const double &mu){
    return log(mu);
  };

  static inline double link_func_inv(const double &eta){
    return std::max(exp(eta), DOUBLE_EPS);
  };

  static inline double variance(const double &mu){
    return mu;
  };

  static inline double d_mu_d_eta(const double &eta){
    return std::max(exp(eta), DOUBLE_EPS);
  };

  static inline double dev_resids(const double &y, const double &mu, const double &w){
    return (y > 0.0) ? w * (y * log(y / mu) - (y - mu)) : mu * w;
  };

  static inline double time_offset(const double &delta_t){
    return log(delta_t);
  };
};

#endif
