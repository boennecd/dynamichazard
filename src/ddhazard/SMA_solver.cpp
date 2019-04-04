#include "../ddhazard.h"
#include "../arma_BLAS_LAPACK.h"
#include "../family.h"

double SMA::compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool is_event, const double length,
    family_base &f){
  double c0 = 0.;
  double c1;

  for(int i = 0; i < 100; i++){
    const double eta = c0 + offset;
    auto trunc_eta = f.truncate_eta(is_event, eta, exp(eta), length);
    const double d1 = f.d_log_like(is_event, trunc_eta, length);
    const double d2 = f.dd_log_like(is_event, trunc_eta, length);

    double &&intermediate =
      (2. * coef1 * c0 + coef2 - w * d1) / (2. * coef1 - w * d2);
    c1 = c0 - intermediate;

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
}

void SMA::solve(){
  double bin_tstop = p_dat.min_start;
  const bool uses_at_risk_length = fam.uses_at_risk_length();
  for (int t = 1; t < p_dat.d + 1; t++){
    const double bin_number = t - 1;
    const double bin_tstart = bin_tstop;
    const double delta_t = p_dat.I_len[t - 1];
    bin_tstop += delta_t;

    // E-step: Prediction step
    p_dat.a_t_less_s.col(t - 1) =
      p_dat.state_trans->map(p_dat.a_t_t_s.col(t - 1)).sv;
    p_dat.V_t_less_s.slice(t - 1) =
      p_dat.state_trans->map(p_dat.V_t_t_s.slice(t - 1)).sv +
      delta_t * p_dat.err_state->map(p_dat.Q).sv;

    if(p_dat.debug){
      std::stringstream str;
      str << t << "|" << t - 1;

      my_print(p_dat, p_dat.a_t_less_s.col(t - 1), "a_(" + str.str() + ")");
      my_print(p_dat, p_dat.V_t_less_s.slice(t - 1), "V_(" + str.str() + ")");
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(p_dat.V_t_less_s.slice(t - 1));
    }

    // E-step: scoring step
    arma::uvec r_set = get_risk_set(p_dat, t);
    arma::vec a(p_dat.a_t_t_s.colptr(t), p_dat.state_dim, false);
    arma::mat V(p_dat.V_t_t_s.slice(t).memptr(), p_dat.state_dim,
                p_dat.state_dim, false);
    a =  p_dat.a_t_less_s.col(t - 1);
    V = p_dat.V_t_less_s.slice(t - 1);

    if(method == "woodbury"){
      for(auto it = r_set.begin(); it != r_set.end(); it++){

        arma::vec x_(p_dat.X.colptr(*it), p_dat.covar_dim, false);
        const double w = p_dat.weights(*it);
        const double offset = p_dat.fixed_effects(*it);

        // TODO: is there a BLAS dsymv for non-square but symetric matrix
        // vector product?
        auto x_in_state_space = p_dat.state_lp_inv->map(x_);
        const arma::vec inter_vec = V * x_in_state_space.sv;

        const double f1 = std::max(
          1./arma::as_scalar(x_in_state_space.sv.t() * inter_vec), 1e-10);
        const double f2 = arma::as_scalar(x_.t() * p_dat.state_lp->map(a).sv);

        const bool is_event = p_dat.is_event_in_bin(*it) == bin_number;

        const double at_risk_length =
          uses_at_risk_length  ?
          get_at_risk_length(
            p_dat.tstop(*it), bin_tstop, p_dat.tstart(*it), bin_tstart) : 0;

        const double c = compute_length(
          offset, f1 / 2., -f2 * f1, w, is_event, at_risk_length, fam);
        double eta = c + offset;
        const double neg_second_d = - w * fam.dd_log_like(
          is_event, eta, exp(eta), at_risk_length);

        a -= (p_dat.LR * (f2 - c) * f1) * inter_vec;
        sym_mat_rank_one_update(
          - neg_second_d  / (1. + neg_second_d / f1), inter_vec, V);
        V = arma::symmatu(V); // TODO: this surely can be done smarter

    }} else if (method == "cholesky"){
      arma::mat L;
      arma::mat L_inv = arma::inv_sympd(V); // only temporary
      symmetric_mat_chol(L_inv, L); // Cholesky decomposition of information matrix
      L_inv = L;
      square_tri_inv(L_inv); // V = L_inv^T * L_inv
      arma::vec inter_vec(L.n_cols);

      for(auto it = r_set.begin(); it != r_set.end(); it++){
        const arma::vec x_(p_dat.X.colptr(*it), p_dat.covar_dim, false);
        const double w = p_dat.weights(*it);
        const double offset = p_dat.fixed_effects(*it);

        tri_mat_times_vec(L_inv, x_, inter_vec, false);

        const double f1 =
          std::max(1./arma::dot(inter_vec, inter_vec), 1e-10);
        const double f2 = arma::dot(x_, p_dat.state_lp->map(a).sv);

        const bool is_event = p_dat.is_event_in_bin(*it) == bin_number;
        const double at_risk_length =
          uses_at_risk_length  ?
          get_at_risk_length(
            p_dat.tstop(*it), bin_tstop, p_dat.tstart(*it), bin_tstart) : 0;

        const double c = compute_length(
          offset, f1 / 2., -f2 * f1, w, is_event, at_risk_length, fam);
        double eta = c + offset;
        const double neg_second_d = - w * fam.dd_log_like(
          is_event, eta, exp(eta), at_risk_length);

        tri_mat_times_vec(L_inv, inter_vec, true);
        a -=  (p_dat.LR * (f2 - c) * f1) * inter_vec;

        arma::vec rank_1_update_vec(x_ * sqrt(neg_second_d));
        rank_1_update_vec = p_dat.state_lp_inv->map(rank_1_update_vec).sv;
        chol_rank_one_update(L, rank_1_update_vec);
        L_inv = L;
        square_tri_inv(L_inv);
      }

      V = L_inv.t() * L_inv;
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
      my_debug_logger(p_dat)
        << "Condition number of V_(" + str.str() + ") is "
        << arma::cond(p_dat.V_t_t_s.slice(t));
    }

    arma::mat V_t_less_s_inv;
    inv_sympd(V_t_less_s_inv, p_dat.V_t_less_s.slice(t - 1), p_dat.use_pinv,
              "ddhazard_fit_cpp estimation error: Failed to invert V_(t|t-1)");
    p_dat.B_s.slice(t - 1) =
      p_dat.state_trans->map(p_dat.V_t_t_s.slice(t - 1), right).sv *
      V_t_less_s_inv;
  }
}
