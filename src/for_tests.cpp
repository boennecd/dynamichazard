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
arma::umat
  sample_n_count_replicas_indices_test(const int size, arma::vec probs)
  {
    std::map<arma::uword, arma::uword> tmp =
      sample_n_count_replicas<sample_indices>(size, probs);
    arma::umat out(tmp.size(), 2L);

    arma::uword j = 0L;
    for (auto it = tmp.begin(); it != tmp.end(); it++, ++j)
    {
      out(j, 0L) = it->first;
      out(j, 1L) = it->second;
    }

    return out;
  }

// [[Rcpp::export]]
arma::umat
  sample_n_count_replicas_systematic_test(const int size, arma::vec probs)
  {
    std::map<arma::uword, arma::uword> tmp =
      sample_n_count_replicas<systematic_resampling>(size, probs);
    arma::umat out(tmp.size(), 2L);

    arma::uword j = 0L;
    for (auto it = tmp.begin(); it != tmp.end(); it++, ++j)
    {
      out(j, 0L) = it->first;
      out(j, 1L) = it->second;
    }

    return out;
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
#include "bigglm_wrapper.h"

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
    logistic fam;
    return(bigglm_updateQR::update(
        qr, X, eta, offset, at_risk_length, y, w, fam));
  } else if (is_exponential_model(model)){
    exponential fam;
    return(bigglm_updateQR::update(
        qr, X, eta, offset, at_risk_length, y, w, fam));
  }
}

// [[Rcpp::export]]
double SMA_hepler_logit_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y){
  logistic f;
  return SMA::compute_length(
    offset, coef1, coef2, w, y, 0., f);
}

// [[Rcpp::export]]
double SMA_hepler_exp_compute_length(
    const double offset, const double coef1, const double coef2,
    const double w, const bool y, const double length){
  exponential f;
  return SMA::compute_length(
    offset, coef1, coef2, w, y, length, f);
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

// [[Rcpp::export]]
arma::mat solve_LU_inv(const arma::mat &A){
  LU_factorization fac(A);
  return fac.solve();
}

// [[Rcpp::export]]
arma::mat solve_LU_mat(const arma::mat &A, const arma::mat &B){
  LU_factorization fac(A);
  return fac.solve(B);
}

// [[Rcpp::export]]
arma::vec solve_LU_vec(const arma::mat &A, const arma::vec &B){
  LU_factorization fac(A);
  return fac.solve(B);
}

// [[Rcpp::export]]
arma::mat qr_qty_mat_test(const arma::mat &A, const arma::mat &B){
  QR_factorization fac(A);
  return fac.qy(B, true);
}

// [[Rcpp::export]]
arma::vec qr_qty_vec_test(const arma::mat &A, const arma::vec &B){
  QR_factorization fac(A);
  return fac.qy(B, true);
}

// [[Rcpp::export]]
arma::mat qr_R_test(const arma::mat &A){
  return QR_factorization(A).R();
}

// [[Rcpp::export]]
arma::mat selection_matrix_map_mat_test(arma::mat L, arma::mat X, bool is_right, bool is_inv){
  selection_matrix S_L(L);
  if(is_inv)
    return S_L.map_inv(X, is_right);

  return S_L.map(X, is_right);
}

// [[Rcpp::export]]
arma::vec selection_matrix_map_vec_test(arma::mat L, arma::vec X, bool is_inv){
  selection_matrix S_L(L);
  if(is_inv)
    return S_L.map_inv(X);

  return S_L.map(X);
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

// -------------------------------------------------- //

#include "lin_maps.h"

// [[Rcpp::export]]
Rcpp::List linear_mapper_test(
    const arma::mat &A, const arma::vec x, const arma::mat X,
    const arma::vec z, const arma::mat Z, std::string type,
    const arma::mat R){
  bool has_R = R.n_elem > 0;
  std::unique_ptr<linear_mapper> ptr;

         if (type == "dens_mapper"){
    ptr.reset(new dens_mapper(A));
  } else if (type == "select_mapper"){
    ptr.reset(new select_mapper(A));
  } else if (type == "inv_mapper"){
    ptr.reset(new inv_mapper(A));
  } else if (type == "inv_sub_mapper"){
    ptr.reset(new inv_sub_mapper(A, R));
  } else
    Rcpp::stop("type not implemented");

  if(!has_R)
    return Rcpp::List::create(
      Rcpp::Named("A") = ptr->map(),

      Rcpp::Named("A_x")     = arma::vec(ptr->map(x       , dont_trans).sv),
      Rcpp::Named("A_T_z")   = arma::vec(ptr->map(z       , trans     ).sv),

      Rcpp::Named("A_X")     = arma::mat(ptr->map(X, left , dont_trans).sv),
      Rcpp::Named("X_A_T")   = arma::mat(ptr->map(X, right, dont_trans).sv),
      Rcpp::Named("A_X_A_T") = arma::mat(ptr->map(X, both , dont_trans).sv),

      Rcpp::Named("A_T_Z")   = arma::mat(ptr->map(Z, left , trans     ).sv),
      Rcpp::Named("Z_A")     = arma::mat(ptr->map(Z, right, trans     ).sv),
      Rcpp::Named("A_T_Z_A") = arma::mat(ptr->map(Z, both , trans     ).sv)
    );

  return Rcpp::List::create(
    Rcpp::Named("A") = ptr->map(),

    Rcpp::Named("A_x")     = arma::vec(ptr->map(x       , dont_trans).sv),

    Rcpp::Named("A_X")     = arma::mat(ptr->map(X, left , dont_trans).sv),
    Rcpp::Named("X_A_T")   = arma::mat(ptr->map(X, right, dont_trans).sv),
    Rcpp::Named("A_X_A_T") = arma::mat(ptr->map(X, both , dont_trans).sv)
  );
}
