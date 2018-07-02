#ifndef BIGGLM_WRAPPER
#define BIGGLM_WRAPPER

#include <memory>
#include "arma_n_rcpp.h"
#include "family.h"

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


class bigglm_updateQR {
  // match logic in update.bigqr
  static arma::vec linkinv(
      const arma::vec&, const arma::vec&, const arma::vec&, family_base&);

  static arma::vec d_mu_d_eta(
      const arma::vec&, const arma::vec&, const arma::vec&, family_base&);

  static arma::vec variance(
      const arma::vec&, const arma::vec&, const arma::vec&, family_base&);

public:
  static void update(
      qr_obj &qr, // Previous/starting value. Will be overwritten
      const arma::mat&, const arma::vec&,
      const arma::vec&, const arma::vec&,
      arma::vec &, const arma::vec&, family_base&);
};

arma::vec bigglm_regcf(qr_obj &qr);

#endif
