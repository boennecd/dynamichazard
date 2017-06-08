#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


// See
//  https://stackoverflow.com/a/42339658/5861244
//  https://ironholds.org/registering-routines/

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .Call calls */
  extern SEXP dynamichazard_bigglm_regcf_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_bigglm_updateQR_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_chol_rank_one_update(SEXP, SEXP);
extern SEXP dynamichazard_ddhazard_fit_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_get_risk_obj_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_IWLS_logit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_IWLS_poisson(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_logLike_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_SMA_hepler_exp_compute_length(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_SMA_hepler_exp_second_d(SEXP, SEXP, SEXP);
extern SEXP dynamichazard_SMA_hepler_logit_compute_length(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dynamichazard_SMA_hepler_logit_second_d(SEXP, SEXP);
extern SEXP dynamichazard_square_tri_inv(SEXP, SEXP);
extern SEXP dynamichazard_sym_mat_rank_one_update(SEXP, SEXP, SEXP);
extern SEXP dynamichazard_symmetric_mat_chol(SEXP, SEXP);
extern SEXP dynamichazard_tri_mat_times_vec(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"dynamichazard_bigglm_regcf_rcpp",               (DL_FUNC) &dynamichazard_bigglm_regcf_rcpp,                6},
  {"dynamichazard_bigglm_updateQR_rcpp",            (DL_FUNC) &dynamichazard_bigglm_updateQR_rcpp,            12},
  {"dynamichazard_chol_rank_one_update",            (DL_FUNC) &dynamichazard_chol_rank_one_update,             2},
  {"dynamichazard_ddhazard_fit_cpp",                (DL_FUNC) &dynamichazard_ddhazard_fit_cpp,                38},
  {"dynamichazard_get_risk_obj_rcpp",               (DL_FUNC) &dynamichazard_get_risk_obj_rcpp,               11},
  {"dynamichazard_IWLS_logit",                      (DL_FUNC) &dynamichazard_IWLS_logit,                       6},
  {"dynamichazard_IWLS_poisson",                    (DL_FUNC) &dynamichazard_IWLS_poisson,                     6},
  {"dynamichazard_logLike_cpp",                     (DL_FUNC) &dynamichazard_logLike_cpp,                     11},
  {"dynamichazard_SMA_hepler_exp_compute_length",   (DL_FUNC) &dynamichazard_SMA_hepler_exp_compute_length,    6},
  {"dynamichazard_SMA_hepler_exp_second_d",         (DL_FUNC) &dynamichazard_SMA_hepler_exp_second_d,          3},
  {"dynamichazard_SMA_hepler_logit_compute_length", (DL_FUNC) &dynamichazard_SMA_hepler_logit_compute_length,  5},
  {"dynamichazard_SMA_hepler_logit_second_d",       (DL_FUNC) &dynamichazard_SMA_hepler_logit_second_d,        2},
  {"dynamichazard_square_tri_inv",                  (DL_FUNC) &dynamichazard_square_tri_inv,                   2},
  {"dynamichazard_sym_mat_rank_one_update",         (DL_FUNC) &dynamichazard_sym_mat_rank_one_update,          3},
  {"dynamichazard_symmetric_mat_chol",              (DL_FUNC) &dynamichazard_symmetric_mat_chol,               2},
  {"dynamichazard_tri_mat_times_vec",               (DL_FUNC) &dynamichazard_tri_mat_times_vec,                4},
  {NULL, NULL, 0}
};

void R_init_dynamichazard(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
