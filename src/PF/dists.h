#ifndef DISTS_H
#define DISTS_H

#include "../arma_n_rcpp.h"
#include "covarmat.h"

/* conditional distribution */
class PF_cdist {
public:
  virtual ~PF_cdist() {}

  /* is it a multivariate normal distribution? */
  virtual bool is_mvn() const = 0;
  /* log density */
  virtual double log_dens(const arma::vec&) const = 0;
  /* gradient of log density */
  virtual arma::vec gradient(const arma::vec&) const = 0;
  /* gradient of log density evaluated at zero vector input */
  virtual arma::vec gradient_zero(const arma::vec&) const = 0;
  /* negative Hessian of log density */
  virtual arma::mat neg_Hessian(const arma::vec&) const = 0;
};

class state_fw : public PF_cdist {
  const arma::vec &parent;
  const arma::mat &F;
  const covarmat &Q;
  const arma::mat QiF;
  const arma::vec mu;

public:
  state_fw(const arma::vec&, const arma::mat&, const covarmat&);
  ~state_fw() = default;

  static double log_dens_func(
      const arma::vec&, const arma::vec&,
      const arma::mat&, const covarmat&);

  bool is_mvn() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec&) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

class state_bw : public PF_cdist {
  const arma::vec &child;
  const arma::mat &F;
  const covarmat &Q;
  const arma::mat FtQi;
  const arma::mat n_hes;

public:
  state_bw(const arma::vec&, const arma::mat&, const covarmat&);
  ~state_bw() = default;

  static double log_dens_func(
      const arma::vec&, const arma::vec&,
      const arma::mat&, const covarmat&);

  bool is_mvn() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec&) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

/*----------------------------------------*/

class artificial_prior : public PF_cdist {
  const arma::vec &mut;
  const covarmat &Qt;
  const arma::vec dz;

public:
  artificial_prior(const arma::vec&, const covarmat&);
  ~artificial_prior() = default;

  bool is_mvn() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec&) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

class artificial_prior_generator {
  using uword = arma::uword;

  const arma::mat &F;
  const covarmat &Q;
  const arma::vec &mu_0;
  const covarmat &Q_0;

  std::map<const uword, const arma::vec> mt;
  std::map<const uword, const covarmat> Pt;

public:
  artificial_prior_generator(const arma::mat&, const covarmat&,
                             const arma::vec&, const covarmat&);

  artificial_prior get_artificial_prior(const arma::uword);
};

/*----------------------------------------*/

class observational_cdist : public PF_cdist {
  const arma::mat &X;
  const arma::vec &y;
  const arma::vec &tstart;
  const arma::vec &tstop;
  const double bin_start;
  const double bin_stop;

public:
  observational_cdist(
    const arma::mat&, const arma::vec&, const arma::vec&,
    const arma::vec&, const double, const double);
  ~observational_cdist() = default;

  bool is_mvn() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec&) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

#endif