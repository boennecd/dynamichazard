// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// chol_rank_one_update
void chol_rank_one_update(arma::mat& R, arma::vec x);
RcppExport SEXP dynamichazard_chol_rank_one_update(SEXP RSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    chol_rank_one_update(R, x);
    return R_NilValue;
END_RCPP
}
// square_tri_inv
void square_tri_inv(const arma::mat& R, arma::mat& out);
RcppExport SEXP dynamichazard_square_tri_inv(SEXP RSEXP, SEXP outSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type out(outSEXP);
    square_tri_inv(R, out);
    return R_NilValue;
END_RCPP
}
// symmetric_mat_chol
void symmetric_mat_chol(const arma::mat& A, arma::mat& out);
RcppExport SEXP dynamichazard_symmetric_mat_chol(SEXP ASEXP, SEXP outSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type out(outSEXP);
    symmetric_mat_chol(A, out);
    return R_NilValue;
END_RCPP
}
// tri_mat_times_vec
void tri_mat_times_vec(arma::mat& A, const arma::vec& x, arma::vec& out, bool is_transpose);
RcppExport SEXP dynamichazard_tri_mat_times_vec(SEXP ASEXP, SEXP xSEXP, SEXP outSEXP, SEXP is_transposeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type out(outSEXP);
    Rcpp::traits::input_parameter< bool >::type is_transpose(is_transposeSEXP);
    tri_mat_times_vec(A, x, out, is_transpose);
    return R_NilValue;
END_RCPP
}
// sym_mat_rank_one_update
void sym_mat_rank_one_update(const double alpha, const arma::vec& x, arma::mat& A);
RcppExport SEXP dynamichazard_sym_mat_rank_one_update(SEXP alphaSEXP, SEXP xSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    sym_mat_rank_one_update(alpha, x, A);
    return R_NilValue;
END_RCPP
}
// bigglm_regcf_rcpp
arma::vec bigglm_regcf_rcpp(arma::vec& D, arma::vec& rbar, arma::vec& thetab, double& ss, bool& checked, arma::vec& tol);
RcppExport SEXP dynamichazard_bigglm_regcf_rcpp(SEXP DSEXP, SEXP rbarSEXP, SEXP thetabSEXP, SEXP ssSEXP, SEXP checkedSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type rbar(rbarSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type thetab(thetabSEXP);
    Rcpp::traits::input_parameter< double& >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< bool& >::type checked(checkedSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(bigglm_regcf_rcpp(D, rbar, thetab, ss, checked, tol));
    return rcpp_result_gen;
END_RCPP
}
// IWLS_logit
arma::vec IWLS_logit(const arma::mat& X, const arma::vec& y, arma::vec beta, const arma::vec offsets, const arma::uword it_max, const double eps);
RcppExport SEXP dynamichazard_IWLS_logit(SEXP XSEXP, SEXP ySEXP, SEXP betaSEXP, SEXP offsetsSEXP, SEXP it_maxSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(IWLS_logit(X, y, beta, offsets, it_max, eps));
    return rcpp_result_gen;
END_RCPP
}
// IWLS_poisson
arma::vec IWLS_poisson(const arma::mat& X, const arma::vec& y, arma::vec beta, const arma::vec offsets, const arma::uword it_max, const double eps);
RcppExport SEXP dynamichazard_IWLS_poisson(SEXP XSEXP, SEXP ySEXP, SEXP betaSEXP, SEXP offsetsSEXP, SEXP it_maxSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(IWLS_poisson(X, y, beta, offsets, it_max, eps));
    return rcpp_result_gen;
END_RCPP
}
// ddhazard_fit_cpp
Rcpp::List ddhazard_fit_cpp(arma::mat& X, arma::mat& fixed_terms, const arma::vec& weights, const arma::vec& tstart, const arma::vec& tstop, const arma::colvec& a_0, const arma::vec& fixed_parems_start, arma::mat Q_0, arma::mat Q, const Rcpp::List& risk_obj, const arma::mat& F_, const double eps_fixed_parems, const int max_it_fixed_params, const arma::uword n_max, const double eps, const arma::uword verbose, const int order_, const bool est_Q_0, const std::string method, Rcpp::Nullable<Rcpp::NumericVector> kappa, Rcpp::Nullable<Rcpp::NumericVector> alpha, Rcpp::Nullable<Rcpp::NumericVector> beta, Rcpp::Nullable<Rcpp::NumericVector> NR_eps, Rcpp::Nullable<Rcpp::NumericVector> LR, const std::string model, const std::string M_step_formulation, const int fixed_effect_chunk_size, const bool debug, const unsigned int NR_it_max, const int n_threads, const double ridge_eps, const int n_fixed_terms_in_state_vec, const bool use_pinv, const std::string criteria, const std::string posterior_version, const signed int GMA_max_rep, const double GMA_NR_eps);
RcppExport SEXP dynamichazard_ddhazard_fit_cpp(SEXP XSEXP, SEXP fixed_termsSEXP, SEXP weightsSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP a_0SEXP, SEXP fixed_parems_startSEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP risk_objSEXP, SEXP F_SEXP, SEXP eps_fixed_paremsSEXP, SEXP max_it_fixed_paramsSEXP, SEXP n_maxSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP order_SEXP, SEXP est_Q_0SEXP, SEXP methodSEXP, SEXP kappaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP NR_epsSEXP, SEXP LRSEXP, SEXP modelSEXP, SEXP M_step_formulationSEXP, SEXP fixed_effect_chunk_sizeSEXP, SEXP debugSEXP, SEXP NR_it_maxSEXP, SEXP n_threadsSEXP, SEXP ridge_epsSEXP, SEXP n_fixed_terms_in_state_vecSEXP, SEXP use_pinvSEXP, SEXP criteriaSEXP, SEXP posterior_versionSEXP, SEXP GMA_max_repSEXP, SEXP GMA_NR_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type fixed_terms(fixed_termsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstart(tstartSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstop(tstopSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fixed_parems_start(fixed_parems_startSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_0(Q_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type risk_obj(risk_objSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F_(F_SEXP);
    Rcpp::traits::input_parameter< const double >::type eps_fixed_parems(eps_fixed_paremsSEXP);
    Rcpp::traits::input_parameter< const int >::type max_it_fixed_params(max_it_fixed_paramsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type n_max(n_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const int >::type order_(order_SEXP);
    Rcpp::traits::input_parameter< const bool >::type est_Q_0(est_Q_0SEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type NR_eps(NR_epsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type LR(LRSEXP);
    Rcpp::traits::input_parameter< const std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const std::string >::type M_step_formulation(M_step_formulationSEXP);
    Rcpp::traits::input_parameter< const int >::type fixed_effect_chunk_size(fixed_effect_chunk_sizeSEXP);
    Rcpp::traits::input_parameter< const bool >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type NR_it_max(NR_it_maxSEXP);
    Rcpp::traits::input_parameter< const int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< const double >::type ridge_eps(ridge_epsSEXP);
    Rcpp::traits::input_parameter< const int >::type n_fixed_terms_in_state_vec(n_fixed_terms_in_state_vecSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_pinv(use_pinvSEXP);
    Rcpp::traits::input_parameter< const std::string >::type criteria(criteriaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type posterior_version(posterior_versionSEXP);
    Rcpp::traits::input_parameter< const signed int >::type GMA_max_rep(GMA_max_repSEXP);
    Rcpp::traits::input_parameter< const double >::type GMA_NR_eps(GMA_NR_epsSEXP);
    rcpp_result_gen = Rcpp::wrap(ddhazard_fit_cpp(X, fixed_terms, weights, tstart, tstop, a_0, fixed_parems_start, Q_0, Q, risk_obj, F_, eps_fixed_parems, max_it_fixed_params, n_max, eps, verbose, order_, est_Q_0, method, kappa, alpha, beta, NR_eps, LR, model, M_step_formulation, fixed_effect_chunk_size, debug, NR_it_max, n_threads, ridge_eps, n_fixed_terms_in_state_vec, use_pinv, criteria, posterior_version, GMA_max_rep, GMA_NR_eps));
    return rcpp_result_gen;
END_RCPP
}
// bigglm_updateQR_rcpp
void bigglm_updateQR_rcpp(arma::vec& D, arma::vec& rbar, arma::vec& thetab, double& ss, bool& checked, arma::vec& tol, std::string model, const arma::mat& X, const arma::vec& eta, const arma::vec& offset, arma::vec& y, const arma::vec& w);
RcppExport SEXP dynamichazard_bigglm_updateQR_rcpp(SEXP DSEXP, SEXP rbarSEXP, SEXP thetabSEXP, SEXP ssSEXP, SEXP checkedSEXP, SEXP tolSEXP, SEXP modelSEXP, SEXP XSEXP, SEXP etaSEXP, SEXP offsetSEXP, SEXP ySEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type rbar(rbarSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type thetab(thetabSEXP);
    Rcpp::traits::input_parameter< double& >::type ss(ssSEXP);
    Rcpp::traits::input_parameter< bool& >::type checked(checkedSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< std::string >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    bigglm_updateQR_rcpp(D, rbar, thetab, ss, checked, tol, model, X, eta, offset, y, w);
    return R_NilValue;
END_RCPP
}
// SMA_hepler_logit_compute_length
double SMA_hepler_logit_compute_length(const double offset, const double coef1, const double coef2, const double w, const bool y);
RcppExport SEXP dynamichazard_SMA_hepler_logit_compute_length(SEXP offsetSEXP, SEXP coef1SEXP, SEXP coef2SEXP, SEXP wSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const double >::type coef1(coef1SEXP);
    Rcpp::traits::input_parameter< const double >::type coef2(coef2SEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const bool >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(SMA_hepler_logit_compute_length(offset, coef1, coef2, w, y));
    return rcpp_result_gen;
END_RCPP
}
// SMA_hepler_logit_second_d
double SMA_hepler_logit_second_d(const double c, const double offset);
RcppExport SEXP dynamichazard_SMA_hepler_logit_second_d(SEXP cSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(SMA_hepler_logit_second_d(c, offset));
    return rcpp_result_gen;
END_RCPP
}
// SMA_hepler_exp_compute_length
double SMA_hepler_exp_compute_length(const double offset, const double coef1, const double coef2, const double w, const bool y, const double length);
RcppExport SEXP dynamichazard_SMA_hepler_exp_compute_length(SEXP offsetSEXP, SEXP coef1SEXP, SEXP coef2SEXP, SEXP wSEXP, SEXP ySEXP, SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const double >::type coef1(coef1SEXP);
    Rcpp::traits::input_parameter< const double >::type coef2(coef2SEXP);
    Rcpp::traits::input_parameter< const double >::type w(wSEXP);
    Rcpp::traits::input_parameter< const bool >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(SMA_hepler_exp_compute_length(offset, coef1, coef2, w, y, length));
    return rcpp_result_gen;
END_RCPP
}
// SMA_hepler_exp_second_d
double SMA_hepler_exp_second_d(const double c, const double offset, const double length);
RcppExport SEXP dynamichazard_SMA_hepler_exp_second_d(SEXP cSEXP, SEXP offsetSEXP, SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< const double >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(SMA_hepler_exp_second_d(c, offset, length));
    return rcpp_result_gen;
END_RCPP
}
// logLike_cpp
std::vector<double> logLike_cpp(const arma::mat& X, const Rcpp::List& risk_obj, const arma::mat& F, const arma::mat& Q_0, arma::mat Q, const arma::mat& a_t_d_s, const arma::vec& tstart, const arma::vec& tstop, const arma::vec& fixed_effects_offsets, const int order_, const std::string model);
RcppExport SEXP dynamichazard_logLike_cpp(SEXP XSEXP, SEXP risk_objSEXP, SEXP FSEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP a_t_d_sSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP fixed_effects_offsetsSEXP, SEXP order_SEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type risk_obj(risk_objSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q_0(Q_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type a_t_d_s(a_t_d_sSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstart(tstartSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstop(tstopSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fixed_effects_offsets(fixed_effects_offsetsSEXP);
    Rcpp::traits::input_parameter< const int >::type order_(order_SEXP);
    Rcpp::traits::input_parameter< const std::string >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(logLike_cpp(X, risk_obj, F, Q_0, Q, a_t_d_s, tstart, tstop, fixed_effects_offsets, order_, model));
    return rcpp_result_gen;
END_RCPP
}
// get_risk_obj_rcpp
Rcpp::List get_risk_obj_rcpp(const Rcpp::NumericVector& start, const Rcpp::NumericVector& stop, const Rcpp::LogicalVector& event, const double& by, const Rcpp::IntegerVector& start_order, const double& max_T, const Rcpp::IntegerVector& order_by_id_and_rev_start, const Rcpp::IntegerVector& id, const double min_start, const Rcpp::NumericVector& event_times_in, const bool& is_for_discrete_model);
RcppExport SEXP dynamichazard_get_risk_obj_rcpp(SEXP startSEXP, SEXP stopSEXP, SEXP eventSEXP, SEXP bySEXP, SEXP start_orderSEXP, SEXP max_TSEXP, SEXP order_by_id_and_rev_startSEXP, SEXP idSEXP, SEXP min_startSEXP, SEXP event_times_inSEXP, SEXP is_for_discrete_modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type stop(stopSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type event(eventSEXP);
    Rcpp::traits::input_parameter< const double& >::type by(bySEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type start_order(start_orderSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_T(max_TSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type order_by_id_and_rev_start(order_by_id_and_rev_startSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type id(idSEXP);
    Rcpp::traits::input_parameter< const double >::type min_start(min_startSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type event_times_in(event_times_inSEXP);
    Rcpp::traits::input_parameter< const bool& >::type is_for_discrete_model(is_for_discrete_modelSEXP);
    rcpp_result_gen = Rcpp::wrap(get_risk_obj_rcpp(start, stop, event, by, start_order, max_T, order_by_id_and_rev_start, id, min_start, event_times_in, is_for_discrete_model));
    return rcpp_result_gen;
END_RCPP
}
