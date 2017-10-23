// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bigglm_regcf_rcpp
arma::vec bigglm_regcf_rcpp(arma::vec& D, arma::vec& rbar, arma::vec& thetab, double& ss, bool& checked, arma::vec& tol);
RcppExport SEXP _dynamichazard_bigglm_regcf_rcpp(SEXP DSEXP, SEXP rbarSEXP, SEXP thetabSEXP, SEXP ssSEXP, SEXP checkedSEXP, SEXP tolSEXP) {
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
// ddhazard_fit_cpp
Rcpp::List ddhazard_fit_cpp(arma::mat& X, arma::mat& fixed_terms, const arma::vec& weights, const arma::vec& tstart, const arma::vec& tstop, const arma::colvec& a_0, const arma::vec& fixed_parems_start, arma::mat Q_0, arma::mat Q, const Rcpp::List& risk_obj, const arma::mat& F_, const double eps_fixed_parems, const int max_it_fixed_params, const arma::uword n_max, const double eps, const arma::uword verbose, const int order_, const bool est_Q_0, const std::string method, Rcpp::Nullable<Rcpp::NumericVector> kappa, Rcpp::Nullable<Rcpp::NumericVector> alpha, Rcpp::Nullable<Rcpp::NumericVector> beta, Rcpp::Nullable<Rcpp::NumericVector> NR_eps, Rcpp::Nullable<Rcpp::NumericVector> LR, const std::string model, const std::string M_step_formulation, const int fixed_effect_chunk_size, const bool debug, const unsigned int NR_it_max, const int n_threads, const double denom_term, const int n_fixed_terms_in_state_vec, const bool use_pinv, const std::string criteria, const std::string posterior_version, const signed int GMA_max_rep, const double GMA_NR_eps, const int EKF_batch_size);
RcppExport SEXP _dynamichazard_ddhazard_fit_cpp(SEXP XSEXP, SEXP fixed_termsSEXP, SEXP weightsSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP a_0SEXP, SEXP fixed_parems_startSEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP risk_objSEXP, SEXP F_SEXP, SEXP eps_fixed_paremsSEXP, SEXP max_it_fixed_paramsSEXP, SEXP n_maxSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP order_SEXP, SEXP est_Q_0SEXP, SEXP methodSEXP, SEXP kappaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP NR_epsSEXP, SEXP LRSEXP, SEXP modelSEXP, SEXP M_step_formulationSEXP, SEXP fixed_effect_chunk_sizeSEXP, SEXP debugSEXP, SEXP NR_it_maxSEXP, SEXP n_threadsSEXP, SEXP denom_termSEXP, SEXP n_fixed_terms_in_state_vecSEXP, SEXP use_pinvSEXP, SEXP criteriaSEXP, SEXP posterior_versionSEXP, SEXP GMA_max_repSEXP, SEXP GMA_NR_epsSEXP, SEXP EKF_batch_sizeSEXP) {
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
    Rcpp::traits::input_parameter< const double >::type denom_term(denom_termSEXP);
    Rcpp::traits::input_parameter< const int >::type n_fixed_terms_in_state_vec(n_fixed_terms_in_state_vecSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_pinv(use_pinvSEXP);
    Rcpp::traits::input_parameter< const std::string >::type criteria(criteriaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type posterior_version(posterior_versionSEXP);
    Rcpp::traits::input_parameter< const signed int >::type GMA_max_rep(GMA_max_repSEXP);
    Rcpp::traits::input_parameter< const double >::type GMA_NR_eps(GMA_NR_epsSEXP);
    Rcpp::traits::input_parameter< const int >::type EKF_batch_size(EKF_batch_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(ddhazard_fit_cpp(X, fixed_terms, weights, tstart, tstop, a_0, fixed_parems_start, Q_0, Q, risk_obj, F_, eps_fixed_parems, max_it_fixed_params, n_max, eps, verbose, order_, est_Q_0, method, kappa, alpha, beta, NR_eps, LR, model, M_step_formulation, fixed_effect_chunk_size, debug, NR_it_max, n_threads, denom_term, n_fixed_terms_in_state_vec, use_pinv, criteria, posterior_version, GMA_max_rep, GMA_NR_eps, EKF_batch_size));
    return rcpp_result_gen;
END_RCPP
}
// sample_indices_test
arma::uvec sample_indices_test(int size, arma::vec probs);
RcppExport SEXP _dynamichazard_sample_indices_test(SEXP sizeSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_indices_test(size, probs));
    return rcpp_result_gen;
END_RCPP
}
// systematic_resampling_test
arma::uvec systematic_resampling_test(const int size, arma::vec probs);
RcppExport SEXP _dynamichazard_systematic_resampling_test(SEXP sizeSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(systematic_resampling_test(size, probs));
    return rcpp_result_gen;
END_RCPP
}
// mvrnorm_test
arma::vec mvrnorm_test(const arma::vec mu, const arma::mat sigma_chol);
RcppExport SEXP _dynamichazard_mvrnorm_test(SEXP muSEXP, SEXP sigma_cholSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type sigma_chol(sigma_cholSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnorm_test(mu, sigma_chol));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_log_test
double dmvnrm_log_test(const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv);
RcppExport SEXP _dynamichazard_dmvnrm_log_test(SEXP xSEXP, SEXP meanSEXP, SEXP sigma_chol_invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type sigma_chol_inv(sigma_chol_invSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_log_test(x, mean, sigma_chol_inv));
    return rcpp_result_gen;
END_RCPP
}
// bigglm_updateQR_rcpp
void bigglm_updateQR_rcpp(arma::vec& D, arma::vec& rbar, arma::vec& thetab, double& ss, bool& checked, arma::vec& tol, std::string model, const arma::mat& X, const arma::vec& eta, const arma::vec& offset, const arma::vec& at_risk_length, arma::vec& y, const arma::vec& w);
RcppExport SEXP _dynamichazard_bigglm_updateQR_rcpp(SEXP DSEXP, SEXP rbarSEXP, SEXP thetabSEXP, SEXP ssSEXP, SEXP checkedSEXP, SEXP tolSEXP, SEXP modelSEXP, SEXP XSEXP, SEXP etaSEXP, SEXP offsetSEXP, SEXP at_risk_lengthSEXP, SEXP ySEXP, SEXP wSEXP) {
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
    Rcpp::traits::input_parameter< const arma::vec& >::type at_risk_length(at_risk_lengthSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    bigglm_updateQR_rcpp(D, rbar, thetab, ss, checked, tol, model, X, eta, offset, at_risk_length, y, w);
    return R_NilValue;
END_RCPP
}
// SMA_hepler_logit_compute_length
double SMA_hepler_logit_compute_length(const double offset, const double coef1, const double coef2, const double w, const bool y);
RcppExport SEXP _dynamichazard_SMA_hepler_logit_compute_length(SEXP offsetSEXP, SEXP coef1SEXP, SEXP coef2SEXP, SEXP wSEXP, SEXP ySEXP) {
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
// SMA_hepler_exp_compute_length
double SMA_hepler_exp_compute_length(const double offset, const double coef1, const double coef2, const double w, const bool y, const double length);
RcppExport SEXP _dynamichazard_SMA_hepler_exp_compute_length(SEXP offsetSEXP, SEXP coef1SEXP, SEXP coef2SEXP, SEXP wSEXP, SEXP ySEXP, SEXP lengthSEXP) {
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
// PF_cloud_to_rcpp_and_back
Rcpp::List PF_cloud_to_rcpp_and_back(const Rcpp::List& rcpp_list);
RcppExport SEXP _dynamichazard_PF_cloud_to_rcpp_and_back(SEXP rcpp_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type rcpp_list(rcpp_listSEXP);
    rcpp_result_gen = Rcpp::wrap(PF_cloud_to_rcpp_and_back(rcpp_list));
    return rcpp_result_gen;
END_RCPP
}
// chol_rank_one_update_test
void chol_rank_one_update_test(arma::mat& R, arma::vec x);
RcppExport SEXP _dynamichazard_chol_rank_one_update_test(SEXP RSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    chol_rank_one_update_test(R, x);
    return R_NilValue;
END_RCPP
}
// square_tri_inv_test
void square_tri_inv_test(const arma::mat& R, arma::mat& out);
RcppExport SEXP _dynamichazard_square_tri_inv_test(SEXP RSEXP, SEXP outSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type out(outSEXP);
    square_tri_inv_test(R, out);
    return R_NilValue;
END_RCPP
}
// symmetric_mat_chol_test
void symmetric_mat_chol_test(const arma::mat& A, arma::mat& out);
RcppExport SEXP _dynamichazard_symmetric_mat_chol_test(SEXP ASEXP, SEXP outSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type out(outSEXP);
    symmetric_mat_chol_test(A, out);
    return R_NilValue;
END_RCPP
}
// tri_mat_times_vec_test
void tri_mat_times_vec_test(arma::mat& A, const arma::vec& x, arma::vec& out, bool is_transpose);
RcppExport SEXP _dynamichazard_tri_mat_times_vec_test(SEXP ASEXP, SEXP xSEXP, SEXP outSEXP, SEXP is_transposeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type out(outSEXP);
    Rcpp::traits::input_parameter< bool >::type is_transpose(is_transposeSEXP);
    tri_mat_times_vec_test(A, x, out, is_transpose);
    return R_NilValue;
END_RCPP
}
// sym_mat_rank_one_update_test
void sym_mat_rank_one_update_test(const double alpha, const arma::vec& x, arma::mat& A);
RcppExport SEXP _dynamichazard_sym_mat_rank_one_update_test(SEXP alphaSEXP, SEXP xSEXP, SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    sym_mat_rank_one_update_test(alpha, x, A);
    return R_NilValue;
END_RCPP
}
// solve_w_precomputed_chol_test
arma::vec solve_w_precomputed_chol_test(const arma::mat& chol_decomp, const arma::vec& B);
RcppExport SEXP _dynamichazard_solve_w_precomputed_chol_test(SEXP chol_decompSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type chol_decomp(chol_decompSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_w_precomputed_chol_test(chol_decomp, B));
    return rcpp_result_gen;
END_RCPP
}
// lambert_W0_test
double lambert_W0_test(const double x);
RcppExport SEXP _dynamichazard_lambert_W0_test(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(lambert_W0_test(x));
    return rcpp_result_gen;
END_RCPP
}
// trunc_eta_exponential_test
Rcpp::List trunc_eta_exponential_test(const double eta, const double at_risk_length, const bool is_event);
RcppExport SEXP _dynamichazard_trunc_eta_exponential_test(SEXP etaSEXP, SEXP at_risk_lengthSEXP, SEXP is_eventSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const double >::type at_risk_length(at_risk_lengthSEXP);
    Rcpp::traits::input_parameter< const bool >::type is_event(is_eventSEXP);
    rcpp_result_gen = Rcpp::wrap(trunc_eta_exponential_test(eta, at_risk_length, is_event));
    return rcpp_result_gen;
END_RCPP
}
// trunc_eta_exponential_test_log_eps
double trunc_eta_exponential_test_log_eps();
RcppExport SEXP _dynamichazard_trunc_eta_exponential_test_log_eps() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(trunc_eta_exponential_test_log_eps());
    return rcpp_result_gen;
END_RCPP
}
// logLike_cpp
std::vector<double> logLike_cpp(const arma::mat& X, const Rcpp::List& risk_obj, const arma::mat& F, const arma::mat& Q_0, arma::mat Q, const arma::mat& a_t_d_s, const arma::vec& tstart, const arma::vec& tstop, const arma::vec& fixed_effects_offsets, const int order_, const std::string model);
RcppExport SEXP _dynamichazard_logLike_cpp(SEXP XSEXP, SEXP risk_objSEXP, SEXP FSEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP a_t_d_sSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP fixed_effects_offsetsSEXP, SEXP order_SEXP, SEXP modelSEXP) {
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
// parallelglm
arma::vec parallelglm(arma::mat& X, /* Not const but will not be touched */     const arma::vec& Ys, std::string family, arma::vec beta0, arma::vec weights, arma::vec offsets, double tol, int nthreads, int it_max, bool trace);
RcppExport SEXP _dynamichazard_parallelglm(SEXP XSEXP, SEXP YsSEXP, SEXP familySEXP, SEXP beta0SEXP, SEXP weightsSEXP, SEXP offsetsSEXP, SEXP tolSEXP, SEXP nthreadsSEXP, SEXP it_maxSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< /* Not const but will not be touched */     const arma::vec& >::type Ys(YsSEXP);
    Rcpp::traits::input_parameter< std::string >::type family(familySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type offsets(offsetsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< int >::type it_max(it_maxSEXP);
    Rcpp::traits::input_parameter< bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelglm(X, Ys, family, beta0, weights, offsets, tol, nthreads, it_max, trace));
    return rcpp_result_gen;
END_RCPP
}
// PF_smooth
Rcpp::List PF_smooth(const int n_fixed_terms_in_state_vec, arma::mat& X, arma::mat& fixed_terms, const arma::vec& tstart, const arma::vec& tstop, const arma::colvec& a_0, arma::mat& Q_0, arma::mat& Q, const arma::mat Q_tilde, const Rcpp::List& risk_obj, const arma::mat& F, const int n_max, const int order, const int n_threads, const int N_fw_n_bw, const int N_smooth, Rcpp::Nullable<Rcpp::NumericVector> forward_backward_ESS_threshold, const int debug, const int N_first, const std::string method, const std::string smoother, const std::string model);
RcppExport SEXP _dynamichazard_PF_smooth(SEXP n_fixed_terms_in_state_vecSEXP, SEXP XSEXP, SEXP fixed_termsSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP a_0SEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP Q_tildeSEXP, SEXP risk_objSEXP, SEXP FSEXP, SEXP n_maxSEXP, SEXP orderSEXP, SEXP n_threadsSEXP, SEXP N_fw_n_bwSEXP, SEXP N_smoothSEXP, SEXP forward_backward_ESS_thresholdSEXP, SEXP debugSEXP, SEXP N_firstSEXP, SEXP methodSEXP, SEXP smootherSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n_fixed_terms_in_state_vec(n_fixed_terms_in_state_vecSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type fixed_terms(fixed_termsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstart(tstartSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstop(tstopSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Q_0(Q_0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Q_tilde(Q_tildeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type risk_obj(risk_objSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F(FSEXP);
    Rcpp::traits::input_parameter< const int >::type n_max(n_maxSEXP);
    Rcpp::traits::input_parameter< const int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< const int >::type N_fw_n_bw(N_fw_n_bwSEXP);
    Rcpp::traits::input_parameter< const int >::type N_smooth(N_smoothSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type forward_backward_ESS_threshold(forward_backward_ESS_thresholdSEXP);
    Rcpp::traits::input_parameter< const int >::type debug(debugSEXP);
    Rcpp::traits::input_parameter< const int >::type N_first(N_firstSEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const std::string >::type smoother(smootherSEXP);
    Rcpp::traits::input_parameter< const std::string >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(PF_smooth(n_fixed_terms_in_state_vec, X, fixed_terms, tstart, tstop, a_0, Q_0, Q, Q_tilde, risk_obj, F, n_max, order, n_threads, N_fw_n_bw, N_smooth, forward_backward_ESS_threshold, debug, N_first, method, smoother, model));
    return rcpp_result_gen;
END_RCPP
}
// compute_summary_stats
Rcpp::List compute_summary_stats(const Rcpp::List& rcpp_list, unsigned int n_threads, const arma::vec& a_0, const arma::mat& Q, const arma::mat& Q_0);
RcppExport SEXP _dynamichazard_compute_summary_stats(SEXP rcpp_listSEXP, SEXP n_threadsSEXP, SEXP a_0SEXP, SEXP QSEXP, SEXP Q_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type rcpp_list(rcpp_listSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n_threads(n_threadsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q_0(Q_0SEXP);
    rcpp_result_gen = Rcpp::wrap(compute_summary_stats(rcpp_list, n_threads, a_0, Q, Q_0));
    return rcpp_result_gen;
END_RCPP
}
// get_risk_obj_rcpp
Rcpp::List get_risk_obj_rcpp(const Rcpp::NumericVector& start, const Rcpp::NumericVector& stop, const Rcpp::LogicalVector& event, const double& by, const Rcpp::IntegerVector& start_order, const double& max_T, const Rcpp::IntegerVector& order_by_id_and_rev_start, const Rcpp::IntegerVector& id, const double min_start, const Rcpp::NumericVector& event_times_in, const bool& is_for_discrete_model);
RcppExport SEXP _dynamichazard_get_risk_obj_rcpp(SEXP startSEXP, SEXP stopSEXP, SEXP eventSEXP, SEXP bySEXP, SEXP start_orderSEXP, SEXP max_TSEXP, SEXP order_by_id_and_rev_startSEXP, SEXP idSEXP, SEXP min_startSEXP, SEXP event_times_inSEXP, SEXP is_for_discrete_modelSEXP) {
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
// round_if_almost_eq
arma::vec round_if_almost_eq(const arma::vec& x, const arma::uvec& x_ord, const arma::vec& boundaries);
RcppExport SEXP _dynamichazard_round_if_almost_eq(SEXP xSEXP, SEXP x_ordSEXP, SEXP boundariesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type x_ord(x_ordSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type boundaries(boundariesSEXP);
    rcpp_result_gen = Rcpp::wrap(round_if_almost_eq(x, x_ord, boundaries));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dynamichazard_bigglm_regcf_rcpp", (DL_FUNC) &_dynamichazard_bigglm_regcf_rcpp, 6},
    {"_dynamichazard_ddhazard_fit_cpp", (DL_FUNC) &_dynamichazard_ddhazard_fit_cpp, 38},
    {"_dynamichazard_sample_indices_test", (DL_FUNC) &_dynamichazard_sample_indices_test, 2},
    {"_dynamichazard_systematic_resampling_test", (DL_FUNC) &_dynamichazard_systematic_resampling_test, 2},
    {"_dynamichazard_mvrnorm_test", (DL_FUNC) &_dynamichazard_mvrnorm_test, 2},
    {"_dynamichazard_dmvnrm_log_test", (DL_FUNC) &_dynamichazard_dmvnrm_log_test, 3},
    {"_dynamichazard_bigglm_updateQR_rcpp", (DL_FUNC) &_dynamichazard_bigglm_updateQR_rcpp, 13},
    {"_dynamichazard_SMA_hepler_logit_compute_length", (DL_FUNC) &_dynamichazard_SMA_hepler_logit_compute_length, 5},
    {"_dynamichazard_SMA_hepler_exp_compute_length", (DL_FUNC) &_dynamichazard_SMA_hepler_exp_compute_length, 6},
    {"_dynamichazard_PF_cloud_to_rcpp_and_back", (DL_FUNC) &_dynamichazard_PF_cloud_to_rcpp_and_back, 1},
    {"_dynamichazard_chol_rank_one_update_test", (DL_FUNC) &_dynamichazard_chol_rank_one_update_test, 2},
    {"_dynamichazard_square_tri_inv_test", (DL_FUNC) &_dynamichazard_square_tri_inv_test, 2},
    {"_dynamichazard_symmetric_mat_chol_test", (DL_FUNC) &_dynamichazard_symmetric_mat_chol_test, 2},
    {"_dynamichazard_tri_mat_times_vec_test", (DL_FUNC) &_dynamichazard_tri_mat_times_vec_test, 4},
    {"_dynamichazard_sym_mat_rank_one_update_test", (DL_FUNC) &_dynamichazard_sym_mat_rank_one_update_test, 3},
    {"_dynamichazard_solve_w_precomputed_chol_test", (DL_FUNC) &_dynamichazard_solve_w_precomputed_chol_test, 2},
    {"_dynamichazard_lambert_W0_test", (DL_FUNC) &_dynamichazard_lambert_W0_test, 1},
    {"_dynamichazard_trunc_eta_exponential_test", (DL_FUNC) &_dynamichazard_trunc_eta_exponential_test, 3},
    {"_dynamichazard_trunc_eta_exponential_test_log_eps", (DL_FUNC) &_dynamichazard_trunc_eta_exponential_test_log_eps, 0},
    {"_dynamichazard_logLike_cpp", (DL_FUNC) &_dynamichazard_logLike_cpp, 11},
    {"_dynamichazard_parallelglm", (DL_FUNC) &_dynamichazard_parallelglm, 10},
    {"_dynamichazard_PF_smooth", (DL_FUNC) &_dynamichazard_PF_smooth, 22},
    {"_dynamichazard_compute_summary_stats", (DL_FUNC) &_dynamichazard_compute_summary_stats, 5},
    {"_dynamichazard_get_risk_obj_rcpp", (DL_FUNC) &_dynamichazard_get_risk_obj_rcpp, 11},
    {"_dynamichazard_round_if_almost_eq", (DL_FUNC) &_dynamichazard_round_if_almost_eq, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_dynamichazard(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
