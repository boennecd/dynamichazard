library(biglm)

if(interactive()){
  library(survival); library(dynamichazard); library(testthat)
  source("C:/Users/boennecd/Dropbox/skole_backup/phd/dynamichazard/R/test_utils.R")
}

sims <- test_sim_func_logit(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = .5, re_draw = T, beta_start = 0,
                            intercept_start = -3, sds = c(.1, rep(1, 3)))



test_that("Only fixed terms yields the same as bigglm", expect_true(F))
test_that("Mixtures versus previous fit for both types of models", expect_true(F))
test_that("Implement UKF", expect_true(F))
test_that("Implement EKF", expect_true(F))
test_that("Test new control variables", expect_true(F))
