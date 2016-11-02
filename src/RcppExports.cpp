// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// bigglm_updateQR_rcpp
void bigglm_updateQR_rcpp(arma::vec& D, arma::vec& rbar, arma::vec& thetab, double& ss, bool& checked, arma::vec& tol, std::string model, const arma::mat& X, const arma::vec& eta, const arma::vec& offset, arma::vec& y);
RcppExport SEXP dynamichazard_bigglm_updateQR_rcpp(SEXP DSEXP, SEXP rbarSEXP, SEXP thetabSEXP, SEXP ssSEXP, SEXP checkedSEXP, SEXP tolSEXP, SEXP modelSEXP, SEXP XSEXP, SEXP etaSEXP, SEXP offsetSEXP, SEXP ySEXP) {
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
    bigglm_updateQR_rcpp(D, rbar, thetab, ss, checked, tol, model, X, eta, offset, y);
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
Rcpp::List ddhazard_fit_cpp(arma::mat& X, arma::mat& fixed_terms, const arma::vec& tstart, const arma::vec& tstop, const arma::colvec& a_0, const arma::vec& fixed_parems_start, arma::mat Q_0, arma::mat Q, const Rcpp::List& risk_obj, const arma::mat& F_, const double eps_fixed_parems, const int max_it_fixed_parems, const arma::uword n_max, const double eps, const arma::uword verbose, const int order_, const bool est_Q_0, const std::string method, Rcpp::Nullable<Rcpp::NumericVector> kappa, Rcpp::Nullable<Rcpp::NumericVector> alpha, Rcpp::Nullable<Rcpp::NumericVector> beta, Rcpp::Nullable<Rcpp::NumericVector> NR_eps, Rcpp::Nullable<Rcpp::NumericVector> LR, const std::string model, const std::string M_step_formulation, const int fixed_effect_chunk_size, const bool debug, const unsigned int NR_it_max);
RcppExport SEXP dynamichazard_ddhazard_fit_cpp(SEXP XSEXP, SEXP fixed_termsSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP a_0SEXP, SEXP fixed_parems_startSEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP risk_objSEXP, SEXP F_SEXP, SEXP eps_fixed_paremsSEXP, SEXP max_it_fixed_paremsSEXP, SEXP n_maxSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP order_SEXP, SEXP est_Q_0SEXP, SEXP methodSEXP, SEXP kappaSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP NR_epsSEXP, SEXP LRSEXP, SEXP modelSEXP, SEXP M_step_formulationSEXP, SEXP fixed_effect_chunk_sizeSEXP, SEXP debugSEXP, SEXP NR_it_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type fixed_terms(fixed_termsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstart(tstartSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstop(tstopSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type fixed_parems_start(fixed_parems_startSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_0(Q_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type risk_obj(risk_objSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F_(F_SEXP);
    Rcpp::traits::input_parameter< const double >::type eps_fixed_parems(eps_fixed_paremsSEXP);
    Rcpp::traits::input_parameter< const int >::type max_it_fixed_parems(max_it_fixed_paremsSEXP);
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
    rcpp_result_gen = Rcpp::wrap(ddhazard_fit_cpp(X, fixed_terms, tstart, tstop, a_0, fixed_parems_start, Q_0, Q, risk_obj, F_, eps_fixed_parems, max_it_fixed_parems, n_max, eps, verbose, order_, est_Q_0, method, kappa, alpha, beta, NR_eps, LR, model, M_step_formulation, fixed_effect_chunk_size, debug, NR_it_max));
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
List get_risk_obj_rcpp(const NumericVector& start, const NumericVector& stop, const LogicalVector& event, const double& by, const IntegerVector& start_order, const double& max_T, const IntegerVector& order_by_id_and_rev_start, const IntegerVector& id, const bool& is_for_discrete_model);
RcppExport SEXP dynamichazard_get_risk_obj_rcpp(SEXP startSEXP, SEXP stopSEXP, SEXP eventSEXP, SEXP bySEXP, SEXP start_orderSEXP, SEXP max_TSEXP, SEXP order_by_id_and_rev_startSEXP, SEXP idSEXP, SEXP is_for_discrete_modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type stop(stopSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type event(eventSEXP);
    Rcpp::traits::input_parameter< const double& >::type by(bySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type start_order(start_orderSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_T(max_TSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type order_by_id_and_rev_start(order_by_id_and_rev_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type id(idSEXP);
    Rcpp::traits::input_parameter< const bool& >::type is_for_discrete_model(is_for_discrete_modelSEXP);
    rcpp_result_gen = Rcpp::wrap(get_risk_obj_rcpp(start, stop, event, by, start_order, max_T, order_by_id_and_rev_start, id, is_for_discrete_model));
    return rcpp_result_gen;
END_RCPP
}
