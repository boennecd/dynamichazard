#ifndef DISTS_H
#define DISTS_H

#include "../arma_n_rcpp.h"
#include "covarmat.h"
#include <memory>

/* conditional distribution */
class PF_cdist {
public:
  virtual ~PF_cdist() = default;

  /* is it a multivariate normal distribution? */
  virtual bool is_mvn() const = 0;
  /* is `gradient_zero` and `neg_Hessian` constant? */
  virtual bool is_grad_z_hes_const() const = 0;
  /* mean of the distribution if the distribution is a distribution in
   * the argument */
  virtual const arma::vec& get_mean() const = 0;
  /* dimension of the present coefficient vector */
  virtual arma::uword dim() const = 0;
  /* log density */
  virtual double log_dens(const arma::vec&) const = 0;
  /* gradient of log density */
  virtual arma::vec gradient(const arma::vec&) const = 0;
  /* gradient of log density evaluated at zero vector input */
  virtual arma::vec gradient_zero(const arma::vec*) const = 0;
  /* negative Hessian of log density */
  virtual arma::mat neg_Hessian(const arma::vec&) const = 0;
};

class state_fw final : public PF_cdist {
  const arma::vec &parent;
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
  bool is_grad_z_hes_const() const override;
  const arma::vec& get_mean() const override;
  arma::uword dim() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec*) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

class state_bw final : public PF_cdist {
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
  bool is_grad_z_hes_const() const override;
  const arma::vec& get_mean() const override;
  arma::uword dim() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec*) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

/*----------------------------------------*/

class artificial_prior final : public PF_cdist {
  const arma::vec &mut;
  const covarmat &Qt;
  const arma::vec dz;

public:
  artificial_prior(const arma::vec&, const covarmat&);
  ~artificial_prior() = default;

  bool is_mvn() const override;
  bool is_grad_z_hes_const() const override;
  const arma::vec& get_mean() const override;
  arma::uword dim() const override;
  double log_dens(const arma::vec&) const override;
  arma::vec gradient(const arma::vec&) const override;
  arma::vec gradient_zero(const arma::vec*) const override;
  arma::mat neg_Hessian(const arma::vec&) const override;
};

class artificial_prior_generator {
  using uword = arma::uword;

  const arma::mat &F;
  const covarmat &Q;

  std::map<const uword, const arma::vec> mt;
  std::map<const uword, const covarmat> Pt;

public:
  artificial_prior_generator(const arma::mat&, const covarmat&,
                             const arma::vec&, const covarmat&);

  artificial_prior get_artificial_prior(const arma::uword);
};

/*----------------------------------------*/

std::shared_ptr<PF_cdist> get_observational_cdist(
    const std::string&, const arma::mat&, const arma::uvec&,
    const arma::vec&, const arma::vec&, const arma::vec&, const double,
    const double,const bool multithreaded = false);

#endif
