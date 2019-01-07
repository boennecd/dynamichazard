#include "ddhazard.h"
#include "family.h"

class logLike_link_term_helper {
protected:
  using uvec_iter = arma::uvec::const_iterator;

  const arma::mat &X;
  const arma::vec &tstart;
  const arma::vec &tstop;
  const arma::ivec &is_event_in_bin;
  const arma::uword n_parems;
  const arma::vec &fixed_effects_offsets;
  family_base &fam;
  const bool uses_at_risk_length;

public:
  logLike_link_term_helper(
    const arma::mat &X_, const arma::vec &tstart_, const arma::vec &tstop_,
    const arma::ivec &is_event_in_bin_,
    const arma::vec &fixed_effects_offsets_, family_base &fam):
    X(X_), tstart(tstart_), tstop(tstop_),
    is_event_in_bin(is_event_in_bin_), n_parems(X.n_rows),
    fixed_effects_offsets(fixed_effects_offsets_), fam(fam),
    uses_at_risk_length(fam.uses_at_risk_length())
  {}

  double link_logLik_terms(
      const arma::vec &a_t,
      const arma::uvec &risk_set,
      const double &bin_start, const double &bin_stop,
      const int bin_number){
    double res(0.0);

    for(uvec_iter it = risk_set.begin(); it != risk_set.end(); ++it){
      const double eta = fixed_effects_offsets(*it) +
        (n_parems > 0 ? arma::dot(a_t.head(n_parems), X.col(*it)) : 0.);
      const double at_risk_length =
        uses_at_risk_length ?
        get_at_risk_length(tstop(*it), bin_stop, tstart(*it), bin_start) :
        0.;

      res += fam.log_like(
        bin_number == is_event_in_bin(*it), eta, exp(eta), at_risk_length);
    }

    return res;
  }
};




// [[Rcpp::export]]
std::vector<double>
  logLike_cpp(const arma::mat &X, const Rcpp::List &risk_obj,
              const arma::mat &F, const arma::mat &Q_0,
              arma::mat Q, const arma::mat &a_t_d_s,
              const arma::vec &tstart, const arma::vec &tstop,
              const arma::vec &fixed_effects_offsets,
              const int order_,
              const std::string model){
  const int d(Rcpp::as<int>(risk_obj["d"]));
  const int n_parems(a_t_d_s.n_rows / order_);
  const bool any_dynamic = n_parems > 0;

  const std::vector<double> &I_len(
      Rcpp::as<std::vector<double> >(risk_obj["I_len"]));

  const arma::ivec &is_event_in_bin(
      Rcpp::as<arma::ivec>(risk_obj["is_event_in"]));

  const Rcpp::List &risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"]));

  double t_0_logLike = 0.;
  arma::mat Q_inv;
  double logLike = 0.;

  if(any_dynamic){
    Q = Q.submat(0, 0, n_parems - 1, n_parems - 1);

    static bool have_fail_to_invert;
    if(!arma::inv_sympd(Q_inv, Q)){
      if(!have_fail_to_invert){
        Rcpp::warning("Failed at least once to invert Q in the log likelihood evaluation with symmetrical inversion method");
        have_fail_to_invert = true;
      }
      if(!arma::inv(Q_inv, Q)){
        Rcpp::stop("Failed to invert Q in log likelihood evalution");
      }
    }

    // compute first terms and normalizing constant
    double log_det_Q, dum;
    arma::log_det(log_det_Q, dum, Q_0);
    t_0_logLike = - 1.0 / 2.0 * log_det_Q
      - n_parems / 2.0 * (log(2.0) + log(M_PI));
  }

  std::unique_ptr<family_base> fam;
  if(model == "logit"){
    fam.reset(new logistic());

  } else if (is_exponential_model(model)){
    fam.reset(new exponential());

  } else if(model == "cloglog"){
    fam.reset(new cloglog());

  } else
    Rcpp::stop("Model '" + model + "' not implemented for logLike method");

  logLike_link_term_helper helper(
      X, tstart, tstop, is_event_in_bin, fixed_effects_offsets,
      *fam.get());

  double log_det_Q, dum;
  arma::log_det(log_det_Q, dum, Q);
  double bin_stop = Rcpp::as<double>(risk_obj["min_start"]);
  for(int t = 1; t <= d; ++t){
    double bin_Start = bin_stop;
    bin_stop += I_len[t - 1];
    arma::vec a_t;

    if(any_dynamic){
      a_t = a_t_d_s.col(t);
      arma::vec delta =
        a_t.head(n_parems) - F.head_rows(n_parems) * a_t_d_s.col(t - 1);

      logLike -=  1.0/ 2.0 * (log_det_Q + log(bin_stop - bin_Start))
        + n_parems / 2.0 * (log(2.0) + log(M_PI));

      logLike -= 1.0 / 2.0 * pow(bin_stop - bin_Start, -1) * arma::as_scalar(
        delta.t() * (Q_inv  * delta));
    }

    const arma::uvec &risk_set = Rcpp::as<arma::uvec>(risk_sets[t - 1]) - 1;

    logLike +=
      helper.link_logLik_terms(a_t, risk_set, bin_Start, bin_stop, t - 1);
  }

  return std::vector<double>{logLike, t_0_logLike};
}
