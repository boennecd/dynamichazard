#ifndef DDFAMILY
#define DDFAMILY

#include "utils.h"
#include <cmath>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define BINOMIAL "binomial"
#define CLOGLOG "cloglog"
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

  virtual double initialize(double) const = 0;

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

  using family_base::log_like;
  double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    double p = linkinv(eta, exp_eta, at_risk_length);
    return outcome ? log(p) : log1p(-p);
  }

  using family_base::d_log_like;
  double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return (exp_eta * (outcome - 1) + outcome) / (exp_eta + 1);
  }

  using family_base::dd_log_like;
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

  double initialize(double y) const override {
    return glm_initialize(y, 1);
  }
};



class cloglog : public virtual family_base, public virtual glm_base {
  const double mu_lb = std::numeric_limits<double>::epsilon();
  const double exp_eta_lb = -log1p(-mu_lb);
  const double eta_lb = log(exp_eta_lb);

  const double mu_ub = 1. - std::numeric_limits<double>::epsilon();
  const double exp_eta_ub = -log(1. - mu_ub);
  const double eta_ub = log(exp_eta_ub);

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
        if(exp_eta < exp_eta_lb){
          ans.eta_trunc = eta_lb;
          ans.exp_eta_trunc = exp_eta_lb;

        } else if (exp_eta > exp_eta_ub){
          ans.eta_trunc = eta_ub;
          ans.exp_eta_trunc = exp_eta_ub;

        } else {
          ans.eta_trunc = eta;
          ans.exp_eta_trunc = exp_eta;

        }

        return ans;
      }

  double linkinv(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
        return -expm1(-exp_eta);
      }

  double mu_eta(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
        return exp(eta - exp_eta);
      }

  double var(
      const double eta, const double exp_eta,
      const double at_risk_length) const override {
        double mu = linkinv(eta, exp_eta, at_risk_length);
        return mu * (1 - mu);
      }

  using family_base::log_like;
  double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
        double p = linkinv(eta, exp_eta, at_risk_length);
        return outcome ? log(p) : log1p(-p);
      }

  using family_base::d_log_like;
  double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
        /*
         * y = 1: \exp\eta / (\exp\exp\eta - 1)
         * y = 0: -\exp\eta
         */
        return outcome ? exp_eta / expm1(exp_eta) : - exp_eta;
      }

  using family_base::dd_log_like;
  double dd_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
        /*
         * y = 1: \frac{\exp\eta(\exp\exp\eta - 1 - \exp(\eta + \exp\eta))}{(\exp\exp\eta - 1)^2}
         *      :   = \frac{-(\exp(-\exp\eta) - 1) - \exp\eta}{(\exp(\exp\eta -\eta) - \exp-\eta)(-(\exp(-\exp\eta) -1)))}
         * The numerator is: 1 - \exp(-\exp\eta) - \exp\eta \overset{y=-\exp\eta}{=} 1 -\exp y + y \approx -(\frac{y^2}{2!}+\frac{y^3}{3!} + \frac{y^4}{4!}+\cdots)
         * y = 0: -\exp\eta
         */
        if(outcome){
          const double m_expexp_eta_m1 = - expm1(-exp_eta),
            d1 = exp(exp_eta - eta) - 1 / exp_eta;

          double numerator = m_expexp_eta_m1 - exp_eta;
          if(eta < -8){
            /* avoid catastrophic cancellation for small eta */
            double y = -exp_eta;
            numerator =
              y * y / 2. * (1. + y / 3. * (1. + y /4. * (1. + y / 5.)));
          }

          return numerator / d1 / m_expexp_eta_m1;
        }

        return - exp_eta;
      }

  /* GLM family like functions */
  double glm_dev_resids(double y, double mu, double wt) const override {
    return - 2 * wt * (y * log(mu) + (1 - y) * log(1 - mu));
  }

  double glm_linkfun(double mu) const override {
    return log(-log1p(-mu));
  }

  double glm_linkinv(double eta) const override {
    const double mu = -expm1(-exp(eta));
    return MAX(
      MIN(mu, 1. -  std::numeric_limits<double>::epsilon()),
      std::numeric_limits<double>::epsilon());
  }

  double glm_variance(double mu) const override {
    return mu * (1 - mu);
  }

  double glm_mu_eta(double eta) const override {
    eta = MIN(eta, 700);
    const double x = exp(eta - exp(eta));
    return MAX(x, std::numeric_limits<double>::epsilon());
  }

  double glm_initialize(double y, double weight) const override {
    return glm_linkfun((weight * y + 0.5)/(weight + 1));
  }

  std::string name() const override {
    return CLOGLOG;
  }

  double initialize(double y) const override {
    return glm_initialize(y, 1);
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

  using family_base::log_like;
  double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return outcome * eta - exp_eta * at_risk_length;
  }

  using family_base::d_log_like;
  double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length) const override {
    return outcome - exp_eta * at_risk_length;
  }

  using family_base::dd_log_like;
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

  double initialize(double y) const override {
    return glm_initialize(y, 1);
  }

  std::string name() const override {
    return POISSON;
  }
};

template<class T>
std::unique_ptr<T> get_fam(const std::string);

#undef MIN
#undef MAX
#endif
