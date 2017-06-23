#include "ddhazard.h"
#ifdef _OPENMP
#include "omp.h"
#endif
#include "arma_utils.h"

inline double GMA_hepler_logit::d1(
    const double eta, const bool is_event, const double at_risk_length){
  if(is_event)
    return(1 / (1 + exp(eta)));


  const double e = exp(eta);
  return(- e / (1 + e));
}

inline double GMA_hepler_logit::d2(
  double eta, const double at_risk_length){
  if(eta < -15. || eta > 15.)
    return(0);

  const double e = exp(eta);
  return - e / pow(1. + e, 2);
}



inline double GMA_hepler_exp::d1(
    const double eta, const bool is_event, const double at_risk_length){
  const double e = exp(eta);
  if(is_event){
    return(1. - e * at_risk_length);
  }

  return(- e * at_risk_length);
}

inline double GMA_hepler_exp::d2(
    double eta, const double at_risk_length){
  return -  exp(eta) * at_risk_length;
}

template<class T>
void GMA<T>::solve(){
  bool have_failed_once = false;
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

      my_print(p_dat, p_dat.a_t_less_s.col(t - 1), "a_(" + str.str() + ")");
      my_print(p_dat, p_dat.V_t_less_s.slice(t - 1), "V_(" + str.str() + ")");
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(p_dat.V_t_less_s.slice(t - 1));
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
              "ddhazard_fit_cpp estimation error: Failed to invert covariance matrix after prediction step");
    arma::vec grad_term = V_t_less_inv * (p_dat.LR * p_dat.a_t_less_s.col(t - 1));

    const arma::vec offsets =
      (p_dat.any_fixed_in_M_step) ?
        p_dat.fixed_terms.cols(r_set).t() * p_dat.fixed_parems :
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

    signed int k;
    for(k = 0; k < max_rep; k++){
      arma::mat X_tilde = X_t;
      arma::vec a_old = a;

      arma::vec eta = (a(*p_dat.span_current_cov).t() * X_t).t() + offsets;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(std::min(p_dat.n_threads, (int)std::ceil(r_set.n_elem / 1000.)))
#endif
      for(arma::uword i = 0; i < r_set.n_elem; i++){
        double w_i = w[i];
        h_1d[i] = w_i * T::d1(eta[i], is_event[i], at_risk_lenght[i]);
        double h_2d_neg = sqrt(- w_i * T::d2(eta[i], at_risk_lenght[i]));
        X_tilde.col(i) *= h_2d_neg;
      }

      /* TODO: This can be faster
       *  1) Do as with EKF use the thread_pool
       */

      X_tilde = arma::symmatu(out_mat_prod(X_tilde));

      if(p_dat.debug){
        my_print(p_dat, X_tilde, "X^T(-p'')X");
        my_debug_logger(p_dat) << "Condition number of X^T(-p'')X is " << arma::cond(X_tilde);
      }

      {
        arma::mat tmp = V_t_less_inv;
        tmp(*p_dat.span_current_cov, *p_dat.span_current_cov) += X_tilde;
        inv_sympd(V, tmp, p_dat.use_pinv,
                  "ddhazard_fit_cpp estimation error: Failed to invert Hessian");
      }

      {
        arma::vec tmp = grad_term;
        if(1. - 1e-15 < p_dat.LR && p_dat.LR < 1. + 1e-15){
          tmp(*p_dat.span_current_cov) +=
            X_tilde * a(*p_dat.span_current_cov) + X_t * h_1d;

        } else {
          tmp(*p_dat.span_current_cov) +=
            X_tilde * a(*p_dat.span_current_cov) + X_t * (h_1d * p_dat.LR);

          tmp += V_t_less_inv * ((1 - p_dat.LR) * a);
        }
        a =  V * tmp;
      }

      if(p_dat.debug){
        std::stringstream str;
        str << "^(" << k + 1 << ")";

        my_print(p_dat, a, "a" + str.str());
        my_print(p_dat, V, "V" + str.str());

        my_debug_logger(p_dat) << "Condition number of V is " << arma::cond(V);
      }

      if(arma::norm(a - a_old, 2) / (arma::norm(a_old, 2) + 1e-8) < NR_eps)
        break;
    }
    if(k == max_rep && !have_failed_once){
      have_failed_once = true;
      std::ostringstream warning;
      warning << "Failed once to make correction step within " << max_rep << " iterations"
              << " and a tolerance of relative change in coefficients of " << NR_eps;
      Rcpp::warning(warning.str());
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

      my_debug_logger(p_dat) << "\n\n_____________________________";

      my_print(p_dat, p_dat.a_t_t_s.col(t), "a_(" + str.str() + ")");
      my_print(p_dat, p_dat.V_t_t_s.slice(t), "V_(" + str.str() + ")\n");
      my_debug_logger(p_dat) << "Condition number of V_(" + str.str() + ") is " << arma::cond(p_dat.V_t_t_s.slice(t));
    }

    p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * V_t_less_inv;
  }
}

// Define classes
template class GMA<GMA_hepler_logit>;
template class GMA<GMA_hepler_exp>;
