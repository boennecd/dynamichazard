#include "ddhazard.h"


inline double Posterior_approx_hepler_logit::NR_delta(
      const double offset, const double coef1, const double coef2,
      const double w, const bool is_event, const double c0){
    const double e = exp(c0 + offset);

    if(is_event){
      return (2. * coef1 * c0 + coef2 - w/(1. + e)) /
        (2. * coef1 + w * e / pow(1. + e, 2));
    }

    return (2. * coef1 * c0 + coef2 + w * e/(1. + e)) /
      (2. * coef1 + w * e / pow(1. + e, 2));
};


inline double Posterior_approx_hepler_logit::get_outcome(
  const bool is_event,
  const double tstart, const double bin_tstart,
  const double tstop, const double bin_tstop){
  return is_event;
};

inline double Posterior_approx_hepler_logit::compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const double y){
  double c0 = 0.;
  double c1;
  bool is_event = y == double_one;

  for(int i = 0; i < 100; i++){
    c1 = c0 - NR_delta(offset, coef1, coef2, w, is_event, c0);

    if(std::abs(c1 - c0) < 1e-5){
      return c1;
    }

    c0 = c1;
  }

  static bool have_failed_once;
  if(!have_failed_once){
    have_failed_once = true;
    Rcpp::warning("Newton Rapshon in prediction step failed at least once\n");
  }

  return c1;
};

inline double Posterior_approx_hepler_logit::second_d(
  const double c, const double offset){
    const double e = exp(c + offset);
    return - e / pow(1. + e, 2);
};

// Export for tests
// [[Rcpp::export]]
double Posterior_approx_hepler_logit_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y){
  return Posterior_approx_hepler_logit::compute_length(
    offset, coef1, coef2, w, y);
};

// [[Rcpp::export]]
double Posterior_approx_hepler_logit_second_d(
    const double c, const double offset){
  return Posterior_approx_hepler_logit::second_d(c, offset);
};






template <class T>
void Posterior_approx<T>::solve(){
  double bin_tstop = p_dat.min_start;

  for (int t = 1; t < p_dat.d + 1; t++){
    const double bin_number = t - 1;
    const double bin_tstart = bin_tstop;
    const double delta_t = p_dat.I_len[t - 1];
    bin_tstop += delta_t;

    // E-step: Prediction step
    p_dat.a_t_less_s.col(t - 1) = p_dat.F_ *  p_dat.a_t_t_s.unsafe_col(t - 1);
    p_dat.V_t_less_s.slice(t - 1) = p_dat.F_ * p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ + delta_t * p_dat.Q;

    if(p_dat.debug){
      std::stringstream str;
      str << t << "|" << t - 1;

      my_print(p_dat.a_t_less_s.col(t - 1), "a_(" + str.str() + ")");
      my_print(p_dat.V_t_less_s.slice(t - 1), "V_(" + str.str() + ")");
    }

    // E-step: scoring step
    arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;

    arma::vec a(p_dat.a_t_t_s.colptr(t), p_dat.space_dim_in_arrays, false);
    arma::mat V(p_dat.V_t_t_s.slice(t).memptr(), p_dat.space_dim_in_arrays,
                p_dat.space_dim_in_arrays, false);
    a =  p_dat.a_t_less_s.col(t - 1);
    V = p_dat.V_t_less_s.slice(t - 1);

    for(auto it = r_set.begin(); it != r_set.end(); it++){

      const arma::vec x_(p_dat.X.colptr(*it), p_dat.n_params_state_vec, false);
      const double w = p_dat.weights(*it);

      const double offset = (p_dat.any_fixed_in_M_step) ?
        arma::dot(p_dat.fixed_parems, p_dat.fixed_terms.col(*it)) : 0.;

      const arma::vec inter_vec =
        V(arma::span::all, p_dat.span_current_cov) * x_;

      const double f1 =
        std::max(1./arma::dot(x_, inter_vec(p_dat.span_current_cov)), 1e-10);
      const double f2 = arma::dot(x_, a.head(p_dat.n_params_state_vec));

      const double y = T::get_outcome(
        p_dat.is_event_in_bin(*it) == bin_number,
        p_dat.tstart(*it), bin_tstart, p_dat.tstop(*it), bin_tstop);

      const double c = T::compute_length(offset, f1 / 2., -f2 * f1, w, y);
      const double neg_second_d = - w * T::second_d(c, offset);

      a -= p_dat.LR * ((f2 - c) * f1) * inter_vec;
      V -= inter_vec *  (inter_vec.t() * (neg_second_d  / (1. + neg_second_d / f1)));

      // TODO: clean-up code from debugging
      // Rcpp::Rcout << "c: "<< c
      //             << "\tw: " << w
      //             << "\t f1/2.: " << f1 / 2.
      //             << "\t -f2 * f1: " << -f2 * f1
      //             << "\t f1: " << f1
      //             << "\t f2: " << f2
      //             << "\tneg 2. deriv: " << neg_second_d
      //             << "\ty: "<< y <<  std::endl;

      // my_print(x_, "x_");
      // my_print(a, "a_");
      // my_print(V, "V");
    }

    if(a.has_inf() || a.has_nan()){
      Rcpp::stop("State vector in correction step has nan or inf elements in in bin " +
        std::to_string(t) + ". Try decreasing the learning rate");

    } else if(V.has_inf() || V.has_nan()){
      Rcpp::stop("Covariance matrix in correction step had inf or nan elements in bin " +
        std::to_string(t) + ". Try decreasing the learning rate");

    }

    if(p_dat.debug){
      std::stringstream str;
      str << t << "|" << t;

      Rcpp::Rcout << "\n\n_____________________________" << std::endl;

      my_print(p_dat.a_t_t_s.col(t), "a_(" + str.str() + ")");
      my_print(p_dat.V_t_t_s.slice(t), "V_(" + str.str() + ")\n");
    }

    arma::mat V_t_less_s_inv;
    inv_sympd(V_t_less_s_inv, p_dat.V_t_less_s.slice(t - 1), p_dat.use_pinv,
              "ddhazard_fit_cpp estimation error: Failed to invert V_(t|t-1)");
    p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * V_t_less_s_inv;
  }
};

// Define classes
template class Posterior_approx<Posterior_approx_hepler_logit>;

