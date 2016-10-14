#include "dynamichazard.h"

class logLike_link_term_helper{
private:
  using uvec_iter = arma::uvec::const_iterator;

  // function that takes in a linear predictor and delta time
  virtual double link_term_if_death(const double &, const double &) = 0;
  virtual double link_term_if_control(const double &, const double &) = 0;

  const arma::mat &X;
  const arma::vec &tstart;
  const arma::vec &tstop;
  const arma::ivec &is_event_in_bin;
  const arma::uvec::elem_type n_parems;

public:
  logLike_link_term_helper(const arma::mat &X_,
                           const arma::vec &tstart_,
                           const arma::vec &tstop_,
                           const arma::ivec &is_event_in_bin_):
    X(X_), tstart(tstart_), tstop(tstop_),
    is_event_in_bin(is_event_in_bin_), n_parems(X.n_rows)
  {}

  double link_logLik_terms(const arma::subview_col<double> a_t,
                           const arma::uvec &risk_set,
                           const double &bin_start, const double &bin_stop,
                           const int bin_number){
    double res(0.0);

    for(uvec_iter it = risk_set.begin(); it != risk_set.end(); ++it){
      const double eta = arma::dot(a_t.head(n_parems), X.col(*it));
      const double delta_t = std::min(tstop(*it), bin_stop) - std::max(tstart(*it), bin_start);

      res += (bin_number == is_event_in_bin(*it)) ?
        link_term_if_death(eta, delta_t) : link_term_if_control(eta, delta_t);
    }

    return res;
  }
};



class logLike_link_term_helper_logit : public logLike_link_term_helper{
private:
  inline double link_term_if_death(const double &eta, const double &delta_t){
    return eta - log(1 + exp(eta));
  };

  inline double link_term_if_control(const double &eta, const double &delta_t){
    return - log(1 + exp(eta));
  };

public:
  logLike_link_term_helper_logit(const arma::mat &X_,
                                 const arma::vec &tstart_,
                                 const arma::vec &tstop_,
                                 const arma::ivec &is_event_in_bin_):
    logLike_link_term_helper(X_, tstart_, tstop_, is_event_in_bin_)
  {}
};



class logLike_link_term_helper_cloglog : public logLike_link_term_helper{
private:
  inline double link_term_if_death(const double &eta, const double &delta_t){
    return eta - exp(eta) * delta_t;
  };

  inline double link_term_if_control(const double &eta, const double &delta_t){
    return - exp(eta) * delta_t;
  };

public:
  logLike_link_term_helper_cloglog(const arma::mat &X_,
                                   const arma::vec &tstart_,
                                   const arma::vec &tstop_,
                                   const arma::ivec &is_event_in_bin_):
    logLike_link_term_helper(X_, tstart_, tstop_, is_event_in_bin_)
  {}
};




// [[Rcpp::export]]
std::vector<double>
  logLike_cpp(const arma::mat &X, const Rcpp::List &risk_obj,
              const arma::mat &F, const arma::mat &Q_0,
              arma::mat Q, const arma::mat &a_t_d_s,
              const arma::vec &tstart, const arma::vec &tstop,
              const int order_ = 1,
              const std::string model = "logit"){
  const int d(Rcpp::as<int>(risk_obj["d"]));
  const int n_parems(a_t_d_s.n_rows / order_);

  const std::vector<double> &I_len(
      Rcpp::as<std::vector<double> >(risk_obj["I_len"]));

  const arma::ivec &is_event_in_bin(
      Rcpp::as<arma::ivec>(risk_obj["is_event_in"]));

  const Rcpp::List &risk_sets(Rcpp::as<Rcpp::List>(risk_obj["risk_sets"]));

  // We want to deal with degenered case
  Q = Q.submat(0, 0, n_parems - 1, n_parems - 1);
  const arma::mat Q_inv = arma::inv_sympd(Q);

  // compute first terms and normalizing constant
  double t_0_logLike = - 1.0 / 2.0 * arma::det(Q_0)
    - n_parems / 2.0 * (log(2.0) + log(M_PI));

  double logLike = - d / 2.0 * arma::det(Q)
    - d * n_parems / 2.0 * (log(2.0) + log(M_PI));

  logLike_link_term_helper *helper;
  if(model == "logit"){
    helper = new logLike_link_term_helper_logit(X, tstart, tstop, is_event_in_bin);
  } else if (model == "exponential"){
    helper = new logLike_link_term_helper_cloglog(X, tstart, tstop, is_event_in_bin);
  } else{
    std::stringstream ss;
    ss << "Model '" << model << "' not implemented for logLike method";
    Rcpp::stop(ss.str());
  }


  double bin_stop = Rcpp::as<double>(risk_obj["min_start"]);
  for(int t = 1; t <= d; ++t){
    const arma::subview_col<double> a_t = a_t_d_s.col(t);
    arma::vec delta = a_t.head(n_parems) - F.head_rows(n_parems) * a_t_d_s.col(t - 1);

    double bin_Start = bin_stop;
    bin_stop += I_len[t - 1];

    logLike -= 1.0 / 2.0 * pow(bin_stop - bin_Start, -1) * arma::as_scalar(
      delta.t() * (Q_inv  * delta));

    const arma::uvec &risk_set = Rcpp::as<arma::uvec>(risk_sets[t - 1]) - 1;

    logLike += helper->link_logLik_terms(a_t, risk_set, bin_Start, bin_stop, t - 1);
  }

  return std::vector<double>{logLike, t_0_logLike};
}
