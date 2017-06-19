# Run before test
# See https://github.com/hadley/testthat/issues/544#issuecomment-260053774

if(interactive()){
  library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("R/test_utils.R")

  exp_model_names <- with(environment(ddhazard), exp_model_names)

  # test_LAPACK_BLAS_wrapper.R
  chol_rank_one_update <- with(environment(ddhazard), chol_rank_one_update)
  square_tri_inv <- with(environment(ddhazard), square_tri_inv)
  symmetric_mat_chol <- with(environment(ddhazard), symmetric_mat_chol)
  tri_mat_times_vec <- with(environment(ddhazard), tri_mat_times_vec)
  sym_mat_rank_one_update <- with(environment(ddhazard), sym_mat_rank_one_update)

  # test_SMA.R
  SMA_hepler_logit_compute_length <-
    with(environment(ddhazard), SMA_hepler_logit_compute_length)
  SMA_hepler_logit_second_d <-
    with(environment(ddhazard), SMA_hepler_logit_second_d)
  SMA_hepler_exp_compute_length <-
    with(environment(ddhazard), SMA_hepler_exp_compute_length)
  SMA_hepler_exp_second_d <-
    with(environment(ddhazard), SMA_hepler_exp_second_d)

  # testbigglm_wrapper.R
  bigglm_updateQR_rcpp <- with(environment(ddhazard), bigglm_updateQR_rcpp)
  bigglm_regcf_rcpp <- with(environment(ddhazard), bigglm_regcf_rcpp)

  # testboot_est.R
  get_frac_n_weights <- with(environment(ddhazard), get_frac_n_weights)

  # testdesign_mat_and_risk_obj.R
  get_permu_data_exp <-  with(environment(ddhazard), get_permu_data_exp)
  get_permu_data_rev_exp <- with(environment(ddhazard), get_permu_data_rev_exp)
  get_order_data_exp <-  with(environment(ddhazard), get_order_data_exp)
  get_order_data_rev_exp <- with(environment(ddhazard), get_order_data_rev_exp)
  get_design_matrix <- environment(ddhazard)$get_design_matrix

  # testIWLS.R
  IWLS_logit <- with(environment(ddhazard), IWLS_logit)

  # teststatic_glm.R
  get_design_matrix <- environment(ddhazard)$get_design_matrix
}

library(testthat)
library(dynamichazard)

options(ddhazard_max_threads = 2)
options(warn=1)

head_neck_cancer <- get_head_neck_cancer_data()
pbc2 <- get_pbc2_data()
options(ddhazard_use_speedglm = F)

# testbigglm_wrapper.R
library(biglm)
bigqr.init <- asNamespace("biglm")$bigqr.init

# testloglike.R and testpredict.R
# library(parallel)
