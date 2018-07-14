#ifndef DDFAMILY
#define DDFAMILY

#include "utils.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define BINOMIAL "binomial"
#define POISSON "poisson"

class glm_base {
public:
  virtual double glm_dev_resids(double, double, double) const = 0;

  virtual double glm_linkfun(double) const = 0;

  virtual double glm_linkinv(double) const = 0;

  virtual double glm_variance(double) const = 0;

  virtual double glm_mu_eta(double) const = 0;

  virtual double glm_initialize(double, double) const = 0;

  virtual std::string name() const = 0;

  // create a virtual, default destructor
  virtual ~glm_base() = default;
};

class family_base {
public:
  virtual bool uses_at_risk_length() const = 0;

  virtual double offset(const double) const = 0;

  virtual trunc_eta_res truncate_eta(
      const bool outcome, const double eta, const double exp_eta,
      const double at_risk_length) const = 0;

  virtual double linkinv(
      const double eta, const double exp_eta,
      const double at_risk_length) const = 0;
  double linkinv(
      const trunc_eta_res res, const double at_risk_length) const {
    return linkinv(res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  };

  virtual double mu_eta(
      const double eta, const double exp_eta,
      const double at_risk_length) const = 0;
  double mu_eta(
      const trunc_eta_res res, const double at_risk_length) const {
    return mu_eta(res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  virtual double var(
      const double eta, const double exp_eta,
      const double at_risk_length) const = 0;
  double var(
      const trunc_eta_res res, const double at_risk_length) const {
    return var(res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  virtual double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const = 0;
  double log_like(
      const bool outcome, const trunc_eta_res res,
      const double at_risk_length) const {
    return log_like(
      outcome, res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  virtual double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const = 0;
  double d_log_like(
      const bool outcome, const trunc_eta_res res,
      const double at_risk_length) const {
    return d_log_like(
      outcome, res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  virtual double dd_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const = 0;
  double dd_log_like(
      const bool outcome, const trunc_eta_res res,
      const double at_risk_length) const {
    return dd_log_like(
      outcome, res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  // create a virtual, default destructor
  virtual ~family_base() = default;
};

class logistic : public virtual family_base, public virtual glm_base {
public:
  bool uses_at_risk_length() const override {
    return false;
  }

  double offset(const double at_risk_length) const override {
    return 0;
  }

  trunc_eta_res truncate_eta(
      const bool outcome, const double eta, const double exp_eta,
      const double at_risk_length) const override {
    trunc_eta_res ans;
    ans.eta_trunc = MIN(MAX(eta, -20), 20);

    ans.exp_eta_trunc = (ans.eta_trunc == eta) ? exp_eta : exp(ans.eta_trunc);
    return ans;
  }

  double linkinv(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
    return 1 / (1 + 1 / exp_eta);
  }

  double mu_eta(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
    double denom = 1 + exp_eta;
    return (exp_eta / denom) / denom;
  }

  double var(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
    return mu_eta(eta, exp_eta, at_risk_length);
  }

  double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    double p = linkinv(eta, exp_eta, at_risk_length);
    return outcome ? log(p) : log(1 - p);
  }

  double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return (exp_eta * (outcome - 1) + outcome) / (exp_eta + 1);
  }

  double dd_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return - mu_eta(eta, exp_eta, at_risk_length);
  }

  /* GLM family like functions */
  double glm_dev_resids(double y, double mu, double wt) const override {
    return - 2 * wt * (y * log(mu) + (1 - y) * log(1 - mu));
  }

  double glm_linkfun(double mu) const override {
    return log(mu / (1 - mu));
  }

  double glm_linkinv(double eta) const override {
    return 1 / (1 + exp(-eta));
  }

  double glm_variance(double mu) const override {
    return mu * (1 - mu);
  }

  double glm_mu_eta(double eta) const override {
    double exp_eta = exp(-eta);
    return (eta < -30 || eta > 30) ?
    std::numeric_limits<double>::epsilon() :
      exp_eta /((1 + exp_eta) * (1 + exp_eta));
  }

  double glm_initialize(double y, double weight) const override {
    return glm_linkfun((weight * y + 0.5)/(weight + 1));
  }

  std::string name() const override {
    return BINOMIAL;
  }
};


class exponential : public virtual family_base, public virtual glm_base  {
public:
  bool uses_at_risk_length() const override {
    return true;
  }

  double offset(const double at_risk_length) const override {
    return std::log(at_risk_length);
  }

  trunc_eta_res truncate_eta(
      const bool outcome, const double eta, const double exp_eta,
      const double at_risk_length) const override {
    return trunc_eta_exponential(outcome, eta, exp_eta, at_risk_length);
  }

  double linkinv(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
    return exp_eta * at_risk_length;
  }

  double mu_eta(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
    return linkinv(eta, exp_eta, at_risk_length);
  }

  double var(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
    return linkinv(eta, exp_eta, at_risk_length);
  }

  double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return outcome * eta - exp_eta * at_risk_length;
  }

  double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return outcome - exp_eta * at_risk_length;
  }

  double dd_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return - exp_eta * at_risk_length;
  }

  /* GLM family like functions. Notice Poisson distribution */
  double glm_dev_resids(double y, double mu, double wt) const override {
    if(y > 0)
      return 2 * (wt * (y * log(y/mu) - (y - mu)));

    return 2 * mu * wt;
  }

  double glm_linkfun(double mu) const override {
    return log(mu);
  }

  double glm_linkinv(double eta) const override {
    return std::max(exp(eta), std::numeric_limits<double>::epsilon());
  }

  double glm_variance(double mu) const override {
    return mu;
  }

  double glm_mu_eta(double eta) const override {
    return std::max(exp(eta), std::numeric_limits<double>::epsilon());
  }

  double glm_initialize(double y, double weight) const override {
    return glm_linkfun(y + 0.1);
  }

  std::string name() const override {
    return POISSON;
  }
};

#undef MIN
#undef MAX
#endif
