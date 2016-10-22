library(biglm)

if(interactive()){
  library(survival); library(dynamichazard); library(testthat)
  source("C:/Users/boennecd/Dropbox/skole_backup/phd/dynamichazard/R/test_utils.R")
}

set.seed(548237)
sims <- test_sim_func_logit(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -4, sds = c(.1, rep(1, 3)))
# sum(sims$res$event)

test_that("Only fixed effects yields same results as bigglm with logit model", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  res1 <- ddhazard(form, data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10,
                   control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3))

  tmp_design <- get_survival_case_Weigths_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10)

  res2 <- bigglm(update(form, Y ~ .), data = tmp_design, family = binomial(), chunksize = 1e3)

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects)))
})



test_that("Only fixed effects yields same results as bigglm with exponential model", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
                   control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3))

  tmp_design <- get_survival_case_Weigths_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10)

  res2 <- bigglm(update(form, Y ~ .), data = tmp_design, family = binomial(), chunksize = 1e3)

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects)))
})


test_that("Only fixed terms yields the same as bigglm", expect_true(F))
test_that("Mixtures versus previous fit for both types of models", expect_true(F))
test_that("Implement UKF", expect_true(F))
test_that("Implement EKF", expect_true(F))
test_that("Test new control variables", expect_true(F))

test_that("Plot accounts for fixed effects", expect_true(F))
test_that("Residuals accounts for fixed effects", expect_true(F))
test_that("Predict accounts for fixed effects", expect_true(F))
test_that("LogLik accounts for fixed effects (inc. df)", expect_true(F))
