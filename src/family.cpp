#include "arma_n_rcpp.h"
#include "family.h"

template<class family>
class family_wrapper {
  static const std::string my_name;

public:
  static std::string name(){
    return my_name;
  }

  static double linkinv(const double eta, const double at_risk_length){
    return family::linkinv(eta, exp(eta), at_risk_length);
  }

  static double mu_eta(const double eta, const double at_risk_length){
    return family::mu_eta(eta, exp(eta), at_risk_length);
  }

  static double var(const double eta, const double at_risk_length){
    return family::var(eta, exp(eta), at_risk_length);
  }

  static double log_like(
      const bool outcome, const double eta, const double at_risk_length){
    return family::log_like(outcome, eta, exp(eta), at_risk_length);
  }

  static double d_log_like(
      const bool outcome, const double eta, const double at_risk_length){
    return family::d_log_like(outcome, eta, exp(eta), at_risk_length);
  }

  static double dd_log_like(
      const bool outcome, const double eta, const double at_risk_length){
    return family::dd_log_like(outcome, eta, exp(eta), at_risk_length);
  }
};

template <>
const std::string family_wrapper<exponential>::my_name = "exponential";
RCPP_MODULE(dd_exponential){
  using namespace Rcpp;

  using wrapper = family_wrapper<exponential>;

  function("name", wrapper::name);

  function("linkinv", wrapper::linkinv,
           List::create(_["eta"], _["at_risk_length"]));
  function("mu_eta", wrapper::mu_eta,
           List::create(_["eta"], _["at_risk_length"]));
  function("var", wrapper::var,
           List::create(_["eta"], _["at_risk_length"]));

  function("log_like", wrapper::log_like,
           List::create(_["outcome"], _["eta"], _["at_risk_length"]));
  function("d_log_like", wrapper::d_log_like,
           List::create(_["outcome"], _["eta"], _["at_risk_length"]));
  function("dd_log_like", wrapper::dd_log_like,
           List::create(_["outcome"], _["eta"], _["at_risk_length"]));
}

template <>
const std::string family_wrapper<logistic>::my_name = "logistic";
RCPP_MODULE(dd_logistic){
  using namespace Rcpp;

  using wrapper = family_wrapper<logistic>;

  function("name", wrapper::name);

  function("linkinv", wrapper::linkinv,
           List::create(_["eta"], _["at_risk_length"]));
  function("mu_eta", wrapper::mu_eta,
           List::create(_["eta"], _["at_risk_length"]));
  function("var", wrapper::var,
           List::create(_["eta"], _["at_risk_length"]));

  function("log_like", wrapper::log_like,
           List::create(_["outcome"], _["eta"], _["at_risk_length"]));
  function("d_log_like", wrapper::d_log_like,
           List::create(_["outcome"], _["eta"], _["at_risk_length"]));
  function("dd_log_like", wrapper::dd_log_like,
           List::create(_["outcome"], _["eta"], _["at_risk_length"]));
}
