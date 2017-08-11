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
#include "estimate_fixed_effects_M_step.h"

using bigglm_updateQR_logit   = bigglm_updateQR<logit_fam>;
using bigglm_updateQR_poisson = bigglm_updateQR<poisson_fam>;

// [[Rcpp::export]]
void bigglm_updateQR_rcpp(arma::vec &D, arma::vec &rbar, arma::vec &thetab,
                          double &ss, bool &checked, arma::vec &tol,
                          std::string model,

                          const arma::mat &X, const arma::vec &eta,
                          const arma::vec &offset, arma::vec &y,
                          const arma::vec &w){
  qr_obj qr;
  qr.D = std::shared_ptr<arma::vec>(&D, [](arma::vec*x) -> void { });
  qr.rbar = std::shared_ptr<arma::vec>(&rbar, [](arma::vec*x) -> void { });
  qr.thetab = std::shared_ptr<arma::vec>(&thetab, [](arma::vec*x) -> void { });
  qr.ss = ss;
  qr.checked = checked;
  qr.tol = std::shared_ptr<arma::vec>(&tol, [](arma::vec*x) -> void { });

  if(model == "logit"){
    return(bigglm_updateQR_logit::update(qr, X, eta, offset, y, w));
  } else if (is_exponential_model(model)){
    return(bigglm_updateQR_poisson::update(qr, X, eta, offset, y, w));
  }
}

// -------------------------------------------------- //

#include "BLAS_and_LAPACK/arma_utils.h"

// [[Rcpp::export]]
arma::vec solve_w_precomputed_chol_test(const arma::mat &chol_decomp, const arma::vec& B){
  return(solve_w_precomputed_chol(chol_decomp, B));
}
