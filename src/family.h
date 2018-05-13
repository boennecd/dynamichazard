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
};

template<class T>
class family_base {
public:
  static inline double linkinv(
      const trunc_eta_res res, const double at_risk_length){
    return T::linkinv(res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  static inline double mu_eta(
      const trunc_eta_res res, const double at_risk_length){
    return T::mu_eta(res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  static inline double var(
      const trunc_eta_res res, const double at_risk_length){
    return T::var(res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  static inline double log_like(
      const bool outcome, const trunc_eta_res res,
      const double at_risk_length){
    return T::log_like(
      outcome, res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  static inline double d_log_like(
      const bool outcome, const trunc_eta_res res,
      const double at_risk_length){
    return T::d_log_like(
      outcome, res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }

  static inline double dd_log_like(
      const bool outcome, const trunc_eta_res res,
      const double at_risk_length){
    return T::dd_log_like(
      outcome, res.eta_trunc, res.exp_eta_trunc, at_risk_length);
  }
};

class logistic : public family_base<logistic>, public glm_base {
public:
  const static bool uses_at_risk_length = false;

  static inline trunc_eta_res truncate_eta(
      const bool outcome, const double eta, const double exp_eta,
      const double at_risk_length){
    trunc_eta_res ans;
    ans.eta_trunc = MIN(MAX(eta, -20), 20);

    ans.exp_eta_trunc = (ans.eta_trunc == eta) ? exp_eta : exp(ans.eta_trunc);
    return ans;
  }

  using family_base<logistic>::linkinv;
  static inline double linkinv(
      const double eta, const double exp_eta, const double at_risk_length){
    return 1 / (1 + 1 / exp_eta);
  }

  using family_base<logistic>::mu_eta;
  static inline double mu_eta(
      const double eta, const double exp_eta, const double at_risk_length){
    double denom = 1 + exp_eta;
    return (exp_eta / denom) / denom;
  }

  using family_base<logistic>::var;
  static inline double var(
      const double eta, const double exp_eta, const double at_risk_length){
    return mu_eta(eta, exp_eta, at_risk_length);
  }

  using family_base<logistic>::log_like;
  static inline double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length){
    double p = linkinv(eta, exp_eta, at_risk_length);
    return outcome ? log(p) : log(1 - p);
  }

  using family_base<logistic>::d_log_like;
  static inline double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length){
    return (exp_eta * (outcome - 1) + outcome) / (exp_eta + 1);
  }

  using family_base<logistic>::dd_log_like;
  static inline double dd_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length){
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


class exponential : public family_base<exponential>, public glm_base {
public:
  const static bool uses_at_risk_length = true;

  static inline trunc_eta_res truncate_eta(
      const bool outcome, const double eta, const double exp_eta,
      const double at_risk_length){
    return trunc_eta_exponential(outcome, eta, exp_eta, at_risk_length);
  }

  using family_base<exponential>::linkinv;
  static inline double linkinv(
      const double eta, const double exp_eta, const double at_risk_length){
    return exp_eta * at_risk_length;
  }

  using family_base<exponential>::mu_eta;
  static inline double mu_eta(
      const double eta, const double exp_eta, const double at_risk_length){
    return linkinv(eta, exp_eta, at_risk_length);
  }

  using family_base<exponential>::var;
  static inline double var(
      const double eta, const double exp_eta, const double at_risk_length){
    return linkinv(eta, exp_eta, at_risk_length);
  }

  using family_base<exponential>::log_like;
  static inline double log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length){
    return outcome * eta - exp_eta * at_risk_length;
  }

  using family_base<exponential>::d_log_like;
  static inline double d_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length){
    return outcome - exp_eta * at_risk_length;
  }

  using family_base<exponential>::dd_log_like;
  static inline double dd_log_like(
      const bool outcome, const double eta,
      const double exp_eta, const double at_risk_length){
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
