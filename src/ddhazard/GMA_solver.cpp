#include "../ddhazard.h"
#include "../utils.h"
#include "../arma_BLAS_LAPACK.h"
#include "../family.h"
#ifdef _OPENMP
#include "omp.h"
#endif

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

      my_print(p_dat, p_dat.a_t_less_s.col(t - 1), "a_(" + str.str() + ")");
      my_print(p_dat, p_dat.V_t_less_s.slice(t - 1), "V_(" + str.str() + ")");
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(p_dat.V_t_less_s.slice(t - 1));
    }

    // E-step: Correction step
    const arma::uvec r_set = get_risk_set(p_dat, t);
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
    arma::vec at_risk_length(r_set.n_elem);
    int i = 0;
    if(T::uses_at_risk_length){
      for(auto it = r_set.begin(); it < r_set.end(); it++, i++){
        at_risk_length[i] = get_at_risk_length(
          p_dat.tstop(*it), bin_tstop, p_dat.tstart(*it), bin_tstart);
      }
    } else
      at_risk_length.zeros();

    arma::vec h_1d(r_set.n_elem);

    unsigned int k, q = X_t.n_rows;
    for(k = 0; k < max_rep; k++){
      arma::mat X_cross(q, q, arma::fill::zeros);
      arma::vec a_old = a;

      arma::vec eta = (a(*p_dat.span_current_cov).t() * X_t).t() + offsets;

#ifdef _OPENMP
      int n_threads = std::max(1, std::min(omp_get_max_threads(), (int)r_set.n_elem / 1000 + 1));
#pragma omp parallel num_threads(n_threads) if(n_threads > 1)
{
#endif
      arma::mat my_X_cross(q, q, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
      for(arma::uword i = 0; i < r_set.n_elem; i++){
        double w_i = w[i];
        double at_i = at_risk_length[i];
        bool is_event_i = is_event[i];

        auto trunc_eta = T::truncate_eta(
          is_event_i, eta[i], exp(eta[i]), at_i);
        h_1d[i] = w_i * T::d_log_like(
          is_event_i, trunc_eta, at_i);
        double h_2d_neg = - w_i * T::dd_log_like(
          is_event_i, trunc_eta, at_i);
        sym_mat_rank_one_update(h_2d_neg, X_t.col(i), my_X_cross);
      }

#ifdef _OPENMP
#pragma omp critical(gma_lock)
{
#endif
      X_cross += my_X_cross;

#ifdef _OPENMP
}
}
#endif

      X_cross = arma::symmatu(X_cross);

      if(p_dat.debug){
        my_print(p_dat, X_cross, "X^T(-p'')X");
        my_debug_logger(p_dat) << "Condition number of X^T(-p'')X is " << arma::cond(X_cross);
      }

      {
        arma::mat tmp = V_t_less_inv;
        tmp(*p_dat.span_current_cov, *p_dat.span_current_cov) += X_cross;
        inv_sympd(V, tmp, p_dat.use_pinv,
                  "ddhazard_fit_cpp estimation error: Failed to invert Hessian");
      }

      {
        arma::vec tmp = grad_term;
        if(1. - 1e-15 < p_dat.LR && p_dat.LR < 1. + 1e-15){
          tmp(*p_dat.span_current_cov) +=
            X_cross * a(*p_dat.span_current_cov) + X_t * h_1d;

        } else {
          tmp(*p_dat.span_current_cov) +=
            X_cross * a(*p_dat.span_current_cov) + X_t * (h_1d * p_dat.LR);

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
      Rcpp::stop("ddhazard_fit_cpp estimation error: State vector in correction step has nan or inf elements in in bin " +
        std::to_string(t) + ". Try decreasing the learning rate");

    } else if(V.has_inf() || V.has_nan()){
      Rcpp::stop("ddhazard_fit_cpp estimation error: Covariance matrix in correction step had inf or nan elements in bin " +
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

template class GMA<logistic>;
template class GMA<exponential>;
