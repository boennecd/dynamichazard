#include "ddhazard.h"

inline double GMA_hepler_logit::d1(
    const double eta, const bool is_event, const double at_risk_length){
  const double e = exp(eta);
  if(is_event){
    return(1 / (1 + e));
  }

  return(- e / (1 + e));
}

inline double GMA_hepler_logit::d2(
  double eta, const double at_risk_length){
  const double e = exp(eta);
  return - e / pow(1. + e, 2);
}

template<class T>
void GMA<T>::solve(){
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

    // E-step: Correction step
    const arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;
    arma::vec a(p_dat.a_t_t_s.colptr(t), p_dat.space_dim_in_arrays, false);
    arma::mat V(p_dat.V_t_t_s.slice(t).memptr(), p_dat.space_dim_in_arrays,
                p_dat.space_dim_in_arrays, false);
    a =  p_dat.a_t_less_s.col(t - 1);
    V = p_dat.V_t_less_s.slice(t - 1);

    arma::mat V_t_less_inv;
    inv_sympd(V_t_less_inv, p_dat.V_t_less_s.slice(t - 1), p_dat.use_pinv,
              "Failed to invert covariance matrix after prediction step");
    arma::vec grad_term = V_t_less_inv * p_dat.a_t_less_s.col(t - 1);

    const arma::vec offsets =
      (p_dat.any_fixed_in_M_step) ?
        p_dat.fixed_parems * p_dat.fixed_terms.cols(r_set) :
        arma::vec(r_set.n_elem, arma::fill::zeros);

    const arma::vec w = p_dat.weights(r_set);
    const arma::mat X_t = p_dat.X.cols(r_set);

    const arma::uvec is_event = p_dat.is_event_in_bin(r_set) == bin_number;
    arma::vec at_risk_lenght(r_set.n_elem);
    int i = 0;
    for(auto it = r_set.begin(); it < r_set.end(); it++, i++){
      at_risk_lenght[i] =
        std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);
    }

    arma::vec h_1d(r_set.n_elem);
    arma::vec h_2d_neg(r_set.n_elem);

    arma::mat X_tilde;
    int k;
    for(k = 1; k < 20; k++){
      arma::vec a_old = a;
      arma::vec eta = (a(p_dat.span_current_cov).t() * X_t).t() + offsets;

      for(arma::uword i = 0; i < r_set.n_elem; i++){
        h_1d[i] = T::d1(eta[i], is_event[i], at_risk_lenght[i]);
        h_2d_neg[i] = - T::d2(eta[i], at_risk_lenght[i]) + p_dat.ridge_eps;
      }
      h_1d %= w;
      h_2d_neg %= w;

      X_tilde = (X_t.each_row() % h_2d_neg.t()) * X_t.t();
      inv_sympd(V, X_tilde + V_t_less_inv, p_dat.use_pinv,
                "Failed to invert Hessian");

      a = V * (X_tilde * a  + grad_term + X_t * h_1d);

      if(arma::norm(a - a_old, 2) / (arma::norm(a_old, 2) + 1e-8) < 1e-4)
        break;
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

    p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * V_t_less_inv;
  }
}

// Define classes
template class GMA<GMA_hepler_logit>;
