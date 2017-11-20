/*
  Rcpp does not search for attributes in sub directories which is why some
  functions are in here. See:
    http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2015-March/008473.html
*/

#include "sample_funcs.h"
#include "PF/densities.h"

// [[Rcpp::export]]
arma::uvec sample_indices_test(int size, arma::vec probs){
  return(sample_indices(size, probs));
}

// [[Rcpp::export]]
arma::uvec systematic_resampling_test(const int size, arma::vec probs){
  return systematic_resampling(size, probs);
}

// [[Rcpp::export]]
arma::vec mvrnorm_test(const arma::vec mu, const arma::mat sigma_chol){
  return mvrnorm(mu, sigma_chol);
}

// [[Rcpp::export]]
double dmvnrm_log_test(
    const arma::vec x, const arma::vec mean, const arma::mat sigma_chol_inv){
  return(dmvnrm_log(x, mean, sigma_chol_inv));
}

// -------------------------------------------------- //

#include "ddhazard.h"
#include "family.h"
#include "estimate_fixed_effects_M_step.h"

using bigglm_updateQR_logit   = bigglm_updateQR<logistic>;
using bigglm_updateQR_poisson = bigglm_updateQR<exponential>;

// [[Rcpp::export]]
void bigglm_updateQR_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                          double &ss, bool &checked, arma::vec &tol,
                          std::string model,

                          const arma::mat &X, const arma::vec &eta,
                          const arma::vec &offset,
                          const arma::vec &at_risk_length,
                          arma::vec &y,
                          const arma::vec &w){
  qr_obj qr;
  qr.D = std::shared_ptr<arma::vec>(&D, [](arma::vec*x) -> void { });
  qr.rbar = std::shared_ptr<arma::vec>(&rbar, [](arma::vec*x) -> void { });
  qr.thetab = std::shared_ptr<arma::vec>(&thetab, [](arma::vec*x) -> void { });
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = std::shared_ptr<arma::vec>(&tol, [](arma::vec*x) -> void { });

  if(model == "logit"){
    return(bigglm_updateQR_logit::update(
        qr, X, eta, offset, at_risk_length, y, w));
  } else if (is_exponential_model(model)){
    return(bigglm_updateQR_poisson::update(
        qr, X, eta, offset, at_risk_length, y, w));
  }
}

// [[Rcpp::export]]
double SMA_hepler_logit_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y){
  return SMA<logistic>::compute_length(
    offset, coef1, coef2, w, y, 0.);
}

// [[Rcpp::export]]
double SMA_hepler_exp_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y, const double length){
  return SMA<exponential>::compute_length(
    offset, coef1, coef2, w, y, length);
}

// -------------------------------------------------- //

#include "PF/PF_utils.h"

// [[Rcpp::export]]
Rcpp::List PF_cloud_to_rcpp_and_back(const Rcpp::List &rcpp_list){
  auto cpp_result = get_clouds_from_rcpp_list(rcpp_list);

  return get_rcpp_list_from_cloud(cpp_result);
}

// -------------------------------------------------- //

#include "arma_BLAS_LAPACK.h"

// [[Rcpp::export]]
void chol_rank_one_update_test(arma::mat &R, arma::vec x){
  return chol_rank_one_update(R, x);
}

// [[Rcpp::export]]
void square_tri_inv_test(const arma::mat &R, arma::mat &out){
  return square_tri_inv(R, out);
}

// [[Rcpp::export]]
void symmetric_mat_chol_test(const arma::mat& A, arma::mat &out){
  return symmetric_mat_chol(A, out);
}

// [[Rcpp::export]]
void tri_mat_times_vec_test(arma::mat &A, const arma::vec &x, arma::vec &out, bool is_transpose){
  return tri_mat_times_vec(A, x, out, is_transpose);
}

// [[Rcpp::export]]
void sym_mat_rank_one_update_test(const double alpha, const arma::vec &x, arma::mat &A){
  return sym_mat_rank_one_update(alpha, x, A);
}

// [[Rcpp::export]]
arma::vec solve_w_precomputed_chol_test(const arma::mat &chol_decomp, const arma::vec& B){
  return(solve_w_precomputed_chol(chol_decomp, B));
}

// -------------------------------------------------- //

#include "utils.h"

// [[Rcpp::export]]
double lambert_W0_test(const double x){
  return lambert_W0(x);
}

// [[Rcpp::export]]
Rcpp::List trunc_eta_exponential_test(
  const double eta, const double at_risk_length, const bool is_event)
{
  auto ans = trunc_eta_exponential(is_event, eta, exp(eta), at_risk_length);

  return Rcpp::List::create(
    Rcpp::Named("eta_trunc") = ans.eta_trunc,
    Rcpp::Named("exp_eta_trunc") = ans.exp_eta_trunc);
}

// [[Rcpp::export]]
double trunc_eta_exponential_test_log_eps(){
  return trunc_eta_exponential_log_eps;
}



