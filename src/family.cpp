#include "arma_n_rcpp.h"
#include "family.h"

/* TODO: make functions that takes in arma::vec and avoid `mapply` */

template<class family>
class family_wrapper {
  static const std::string my_name;

  static bool check_eta_n_r_length(
      const arma::vec &eta, const arma::vec &at_risk_length){
    if((eta.n_elem != at_risk_length.n_elem and at_risk_length.n_elem > 1) or
         eta.n_elem < at_risk_length.n_elem)
      Rcpp::stop("Invalid `eta` and `at_risk_length`");

    return true;
  }

  static bool check_outcome_eta_n_r_length(
      const Rcpp::LogicalVector outcome, const arma::vec &eta,
      const arma::vec &at_risk_length)
    {
      if(outcome.size() != eta.n_elem or
           (eta.n_elem != at_risk_length.n_elem and at_risk_length.n_elem > 1) or
           eta.n_elem < at_risk_length.n_elem)
        Rcpp::stop("Invalid `outcome`, `eta` and `at_risk_length`");

      return true;
    }

public:
  static std::string name(){
    return my_name;
  }

  static Rcpp::NumericVector
  linkinv(const arma::vec &eta, const arma::vec &at_risk_length){
    check_eta_n_r_length(eta, at_risk_length);

    family fam;
    Rcpp::NumericVector out(eta.n_elem);
    const double *e = eta.begin(), *a = at_risk_length.begin();
    bool do_inc_a = at_risk_length.n_elem > 1L;
    for(auto o = out.begin(); o != out.end(); ++o, ++e, a += do_inc_a)
      *o = fam.linkinv(*e, exp(*e), *a);

    return out;
  }

  static Rcpp::NumericVector
  mu_eta(const arma::vec &eta, const arma::vec &at_risk_length){
    check_eta_n_r_length(eta, at_risk_length);

    family fam;
    Rcpp::NumericVector out(eta.n_elem);
    const double *e = eta.begin(), *a = at_risk_length.begin();
    bool do_inc_a = at_risk_length.n_elem > 1L;
    for(auto o = out.begin(); o != out.end(); ++o, ++e, a += do_inc_a)
      *o = fam.mu_eta(*e, exp(*e), *a);

    return out;
  }

  static Rcpp::NumericVector
  var(const arma::vec &eta, const arma::vec &at_risk_length){
    check_eta_n_r_length(eta, at_risk_length);

    family fam;
    Rcpp::NumericVector out(eta.n_elem);
    const double *e = eta.begin(), *a = at_risk_length.begin();
    bool do_inc_a = at_risk_length.n_elem > 1L;
    for(auto o = out.begin(); o != out.end(); ++o, ++e, a += do_inc_a)
      *o = fam.var(*e, exp(*e), *a);

    return out;
  }

  static Rcpp::NumericVector
  log_like(const Rcpp::LogicalVector outcome, const arma::vec &eta,
           const arma::vec &at_risk_length){
    family fam;
    Rcpp::NumericVector out(eta.n_elem);
    const double *e = eta.begin(), *a = at_risk_length.begin();
    const int *r = outcome.begin(), do_inc_a = at_risk_length.n_elem > 1L;
    for(auto o = out.begin(); o != out.end(); ++o, ++r, ++e, a += do_inc_a)
      *o = fam.log_like(*r, *e, exp(*e), *a);

    return out;
  }

  static Rcpp::NumericVector
  d_log_like(const Rcpp::LogicalVector outcome, const arma::vec &eta,
             const arma::vec &at_risk_length){
    family fam;
    Rcpp::NumericVector out(eta.n_elem);
    const double *e = eta.begin(), *a = at_risk_length.begin();
    const int *r = outcome.begin(), do_inc_a = at_risk_length.n_elem > 1L;
    for(auto o = out.begin(); o != out.end(); ++o, ++r, ++e, a += do_inc_a)
      *o = fam.d_log_like(*r, *e, exp(*e), *a);

    return out;
  }

  static Rcpp::NumericVector
  dd_log_like(const Rcpp::LogicalVector outcome, const arma::vec &eta,
              const arma::vec &at_risk_length){
    family fam;
    Rcpp::NumericVector out(eta.n_elem);
    const double *e = eta.begin(), *a = at_risk_length.begin();
    const int *r = outcome.begin(), do_inc_a = at_risk_length.n_elem > 1L;
    for(auto o = out.begin(); o != out.end(); ++o, ++r, ++e, a += do_inc_a)
      *o = fam.dd_log_like(*r, *e, exp(*e), *a);

    return out;
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

template <>
const std::string family_wrapper<cloglog>::my_name = "cloglog";
RCPP_MODULE(dd_cloglog){
  using namespace Rcpp;

  using wrapper = family_wrapper<cloglog>;

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



template<class T>
std::unique_ptr<T> get_fam(const std::string family){
  if(family == BINOMIAL){
    return std::unique_ptr<T>(new logistic());

  } else if(family == POISSON){
    return std::unique_ptr<T>(new exponential());

  } else if(family == CLOGLOG){
    return std::unique_ptr<T>(new cloglog());

  } else
    Rcpp::stop("Family not implemented");
}

template std::unique_ptr<glm_base> get_fam<glm_base>(const std::string);
template std::unique_ptr<family_base> get_fam<family_base>(const std::string);
