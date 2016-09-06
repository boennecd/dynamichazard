// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ddhazard_fit_cpp_prelim
Rcpp::List ddhazard_fit_cpp_prelim(const Rcpp::NumericMatrix& X, const arma::vec& tstart, const arma::vec& tstop, const arma::colvec& a_0, arma::mat Q_0, arma::mat Q, const Rcpp::List& risk_obj, const arma::mat& F_, const int n_max, const double eps, const bool verbose, const bool save_all_output, const int order_, const bool est_Q_0, const std::string method, Rcpp::Nullable<Rcpp::NumericVector> k);
RcppExport SEXP dynamichazard_ddhazard_fit_cpp_prelim(SEXP XSEXP, SEXP tstartSEXP, SEXP tstopSEXP, SEXP a_0SEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP risk_objSEXP, SEXP F_SEXP, SEXP n_maxSEXP, SEXP epsSEXP, SEXP verboseSEXP, SEXP save_all_outputSEXP, SEXP order_SEXP, SEXP est_Q_0SEXP, SEXP methodSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstart(tstartSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type tstop(tstopSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q_0(Q_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type risk_obj(risk_objSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F_(F_SEXP);
    Rcpp::traits::input_parameter< const int >::type n_max(n_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type save_all_output(save_all_outputSEXP);
    Rcpp::traits::input_parameter< const int >::type order_(order_SEXP);
    Rcpp::traits::input_parameter< const bool >::type est_Q_0(est_Q_0SEXP);
    Rcpp::traits::input_parameter< const std::string >::type method(methodSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type k(kSEXP);
    __result = Rcpp::wrap(ddhazard_fit_cpp_prelim(X, tstart, tstop, a_0, Q_0, Q, risk_obj, F_, n_max, eps, verbose, save_all_output, order_, est_Q_0, method, k));
    return __result;
END_RCPP
}
// gen_kalman_filter_cpp
List gen_kalman_filter_cpp(const arma::colvec& a_0, const arma::mat& Q_0, const arma::mat& Q, const arma::mat& F_, const List& risk_sets, const NumericVector& I_len, const int d, const NumericMatrix& X, const IntegerVector& start, const arma::ivec& stop, const arma::ivec& is_event_in_bin, const int& order_);
RcppExport SEXP dynamichazard_gen_kalman_filter_cpp(SEXP a_0SEXP, SEXP Q_0SEXP, SEXP QSEXP, SEXP F_SEXP, SEXP risk_setsSEXP, SEXP I_lenSEXP, SEXP dSEXP, SEXP XSEXP, SEXP startSEXP, SEXP stopSEXP, SEXP is_event_in_binSEXP, SEXP order_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::colvec& >::type a_0(a_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q_0(Q_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type F_(F_SEXP);
    Rcpp::traits::input_parameter< const List& >::type risk_sets(risk_setsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type I_len(I_lenSEXP);
    Rcpp::traits::input_parameter< const int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type stop(stopSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type is_event_in_bin(is_event_in_binSEXP);
    Rcpp::traits::input_parameter< const int& >::type order_(order_SEXP);
    __result = Rcpp::wrap(gen_kalman_filter_cpp(a_0, Q_0, Q, F_, risk_sets, I_len, d, X, start, stop, is_event_in_bin, order_));
    return __result;
END_RCPP
}
// get_risk_obj_rcpp
List get_risk_obj_rcpp(const NumericVector& start, const NumericVector& stop, const LogicalVector& event, const double& by, const IntegerVector& start_order, const double& max_T, const IntegerVector& order_by_id_and_rev_start, const IntegerVector& id);
RcppExport SEXP dynamichazard_get_risk_obj_rcpp(SEXP startSEXP, SEXP stopSEXP, SEXP eventSEXP, SEXP bySEXP, SEXP start_orderSEXP, SEXP max_TSEXP, SEXP order_by_id_and_rev_startSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericVector& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type stop(stopSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type event(eventSEXP);
    Rcpp::traits::input_parameter< const double& >::type by(bySEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type start_order(start_orderSEXP);
    Rcpp::traits::input_parameter< const double& >::type max_T(max_TSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type order_by_id_and_rev_start(order_by_id_and_rev_startSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type id(idSEXP);
    __result = Rcpp::wrap(get_risk_obj_rcpp(start, stop, event, by, start_order, max_T, order_by_id_and_rev_start, id));
    return __result;
END_RCPP
}
