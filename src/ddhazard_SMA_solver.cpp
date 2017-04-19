#include "ddhazard.h"
#include "arma_utils.h"

inline double SMA_hepler_logit::NR_delta(
      const double offset, const double coef1, const double coef2,
      const double w, const double c0, const bool is_event,
      const double length){
    const double e = exp(c0 + offset);

    if(is_event){
      return (2. * coef1 * c0 + coef2 - w/(1. + e)) /
        (2. * coef1 + w * e / pow(1. + e, 2));
    }

    return (2. * coef1 * c0 + coef2 + w * e/(1. + e)) /
      (2. * coef1 + w * e / pow(1. + e, 2));
};

inline double SMA_hepler_logit::compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool is_event, const double length){
  double c0 = 0.;
  double c1;

  for(int i = 0; i < 100; i++){
    c1 = c0 - NR_delta(offset, coef1, coef2, w, c0, is_event, 0.);

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

inline double SMA_hepler_logit::second_d(
  const double c, const double offset, const double length){
    const double e = exp(c + offset);
    return - e / pow(1. + e, 2);
};

// Export for tests
// [[Rcpp::export]]
double SMA_hepler_logit_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y){
  return SMA_hepler_logit::compute_length(
    offset, coef1, coef2, w, y, 0.);
};

// [[Rcpp::export]]
double SMA_hepler_logit_second_d(
    const double c, const double offset){
  return SMA_hepler_logit::second_d(c, offset,0.);
};






inline double SMA_hepler_exp::NR_delta(
    const double offset, const double coef1, const double coef2,
    const double w, const double c0, const bool is_event,
    const double length){
  const double e = exp(c0 + offset + log(length));

  if(is_event){
    return (2. * coef1 * c0 + coef2 - w * (1. - e)) / (2. * coef1 + w * e);
  }

  return (2. * coef1 * c0 + coef2 + w *  e) / (2. * coef1 + w * e);
};

inline double SMA_hepler_exp::compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool is_event, const double length){
  double c0 = 0.;
  double c1;

  for(int i = 0; i < 100; i++){
    c1 = c0 - NR_delta(offset, coef1, coef2, w, c0, is_event, length);

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

inline double SMA_hepler_exp::second_d(
    const double c, const double offset, const double length){
  return -  exp(c + offset + log(length));
};

// Export for tests
// [[Rcpp::export]]
double SMA_hepler_exp_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y, const double length){
  return SMA_hepler_exp::compute_length(
    offset, coef1, coef2, w, y, length);
};

// [[Rcpp::export]]
double SMA_hepler_exp_second_d(
    const double c, const double offset, const double length){
  return SMA_hepler_exp::second_d(c, offset, length);
};








template <class T>
void SMA<T>::solve(){
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
      Rcpp::Rcout << "Condition number of V_(" + str.str() + ") is "
                  << arma::cond(p_dat.V_t_less_s.slice(t - 1)) << std::endl;
    }

    // E-step: scoring step
#ifdef USE_OPEN_BLAS
    // adverse effect on performance unless the dimension of state vector is
    // quite large
    openblas_set_num_threads(1);
#endif

    arma::uvec r_set = Rcpp::as<arma::uvec>(p_dat.risk_sets[t - 1]) - 1;
    arma::vec a(p_dat.a_t_t_s.colptr(t), p_dat.space_dim_in_arrays, false);
    arma::mat V(p_dat.V_t_t_s.slice(t).memptr(), p_dat.space_dim_in_arrays,
                p_dat.space_dim_in_arrays, false);
    a =  p_dat.a_t_less_s.col(t - 1);
    V = p_dat.V_t_less_s.slice(t - 1);

    if(method == "woodbury"){
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

        const bool is_event = p_dat.is_event_in_bin(*it) == bin_number;
        const double at_risk_lenght =
          std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);

        const double c = T::compute_length(offset, f1 / 2., -f2 * f1, w, is_event, at_risk_lenght);
        const double neg_second_d = - w * T::second_d(c, offset, at_risk_lenght);

        a -= p_dat.LR * ((f2 - c) * f1) * inter_vec;
        sym_mat_rank_one_update(
          - neg_second_d  / (1. + neg_second_d / f1), inter_vec, V);

    }} else if (method == "cholesky"){
      arma::mat L;
      arma::mat L_inv = arma::inv_sympd(V); // only temporary
      symmetric_mat_chol(L_inv, L); // Cholesky decomposition of information matrix
      square_tri_inv(L, L_inv); // V = L_inv^T * L_inv
      arma::vec inter_vec(L.n_cols);

      for(auto it = r_set.begin(); it != r_set.end(); it++){
        const arma::vec x_(p_dat.X.colptr(*it), p_dat.n_params_state_vec, false);
        const double w = p_dat.weights(*it);

        const double offset = (p_dat.any_fixed_in_M_step) ?
          arma::dot(p_dat.fixed_parems, p_dat.fixed_terms.col(*it)) : 0.;

        tri_mat_times_vec(L_inv, x_, inter_vec, false);

        const double f1 =
          std::max(1./arma::dot(inter_vec, inter_vec), 1e-10);
        const double f2 = arma::dot(x_, a.head(p_dat.n_params_state_vec));

        const bool is_event = p_dat.is_event_in_bin(*it) == bin_number;
        const double at_risk_lenght =
          std::min(p_dat.tstop(*it), bin_tstop) - std::max(p_dat.tstart(*it), bin_tstart);

        const double c = T::compute_length(offset, f1 / 2., -f2 * f1, w, is_event, at_risk_lenght);
        const double neg_second_d = - w * T::second_d(c, offset, at_risk_lenght);

        tri_mat_times_vec(L_inv, inter_vec, true);
        a -=  (p_dat.LR * (f2 - c) * f1) * inter_vec;

        arma::vec rank_1_update_vec(p_dat.space_dim_in_arrays, arma::fill::zeros);
        rank_1_update_vec(p_dat.span_current_cov) = x_ * sqrt(neg_second_d);
        chol_rank_one_update(L, rank_1_update_vec);
        square_tri_inv(L, L_inv);
      }

      V = L_inv.t() * L_inv;
    }

#ifdef USE_OPEN_BLAS
    openblas_set_num_threads(p_dat.n_threads);
#endif

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
      Rcpp::Rcout << "Condition number of V_(" + str.str() + ") is "
                  << arma::cond(p_dat.V_t_t_s.slice(t)) << std::endl;
    }

    arma::mat V_t_less_s_inv;
    inv_sympd(V_t_less_s_inv, p_dat.V_t_less_s.slice(t - 1), p_dat.use_pinv,
              "ddhazard_fit_cpp estimation error: Failed to invert V_(t|t-1)");
    p_dat.B_s.slice(t - 1) = p_dat.V_t_t_s.slice(t - 1) * p_dat.T_F_ * V_t_less_s_inv;
  }
};

// Define classes
template class SMA<SMA_hepler_logit>;
template class SMA<SMA_hepler_exp>;

