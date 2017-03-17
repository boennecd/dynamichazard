if(interactive()){
  library(dynamichazard); library(testthat); library(survival)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")

  exp_model_names <- with(environment(ddhazard), exp_model_names)
}


# Had issues with win builder. Thus, these lines
test_name <- "UKF"
cat("\nRunning", test_name, "\n")

test_that("UKF throws error when first sigm points weight is zero", {
  expect_error({
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(1, 2), a_0 = c(-3, 0),
      Q = diag(1e-1, 2),
      control = list(kappa = 0, method = "UKF"),
      max_T = 30,
      id = head_neck_cancer$id, order = 1)
  }, regexp = "UKF not implemented for hyperparameters that yield zero weight on first sigma point")
})


test_that("UKF on head_neck works with logit model", {
  suppressMessages(result <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, Q_0 = diag(1, 2), a_0 = c(-3, 0),
    Q = diag(1e-1, 2),
    control = list(est_Q_0 = F, n_max = 10^4, eps = 10^-3,
                   method = "UKF", save_data = F, save_risk_set = F,
                   beta = 0, alpha = 1),
    max_T = 30,
    id = head_neck_cancer$id, order = 1,
    verbose = F,
    model = "logit"
  ))

  # plot(result)
  # save_to_test(result, "UKF1")

  expect_equal(result, read_to_test("UKF1"))
})


set.seed(2972)
sims <- test_sim_func_logit(n_series = 5e2, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = -.5, re_draw = T, beta_start = 1,
                            intercept_start = -1, sds = c(.1, rep(1, 3)))

test_that("UKF does not fail and both methods give the same",{
  res_new <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                  by = 1,
                  data = sims$res,
                  a_0 = rep(0, ncol(sims$res) + 1 - 4),
                  Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                  verbose = F,
                  control = list(eps = 1e-2, est_Q_0 = F, method = "UKF",
                                 kappa = 0.1, alpha = 1, beta = 0),
                  id = sims$res$id,
                  max_T = 10)

  res_old <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                      by = 1,
                      data = sims$res,
                      a_0 = rep(0, ncol(sims$res) + 1 - 4),
                      Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                      control = list(est_Q_0 = F, kappa = 0.1, alpha = 1, beta = 0, eps = 1e-2, method = "UKF_org"),
                      id = sims$res$id,
                      max_T = 10)

  expect_equal(res_new$state_vecs, res_old$state_vecs)
  expect_equal(res_new$state_vars, res_old$state_vars)
})

test_that("Chaning time scale in UKF does no change results when other parems are changed accoridngly",{
  arg_list <- list(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                   by = 1,
                   data = sims$res,
                   a_0 = rep(0, ncol(sims$res) + 1 - 4),
                   Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                   Q = diag(rep(1e-2, ncol(sims$res) + 1 - 4)),
                   control = list(kappa = 0.1, alpha = 1, beta = 0,
                                  est_Q_0 = F, method = "UKF", eps = 1e-2),
                   id = sims$res$id,
                   verbose = F,
                   max_T = 10)

  res <- do.call(ddhazard, arg_list)

  sim_tmp <- sims
  t_mult <- 2.5
  sim_tmp$res$tstart <- sim_tmp$res$tstart * t_mult
  sim_tmp$res$tstop <- sim_tmp$res$tstop * t_mult
  arg_list$data <- sim_tmp$res
  arg_list$by <- t_mult
  arg_list$max_T <- arg_list$max_T * t_mult
  arg_list$Q <- arg_list$Q / t_mult

  res_new_time <- do.call(ddhazard, arg_list)

  expect_equal(res$state_vecs, res_new_time$state_vecs)
  expect_equal(res$state_vars, res_new_time$state_vars)
  expect_equal(res$Q, res_new_time$Q * t_mult)
})

test_that("Testing UKF against prev computed values",{
  arg_list <- list(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                   by = 1,
                   data = sims$res,
                   a_0 = rep(0, ncol(sims$res) + 1 - 4),
                   Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                   Q = diag(rep(1e-2, ncol(sims$res) + 1 - 4)),
                   control = list(kappa = 0.1, alpha = 1, beta = 0,
                                  est_Q_0 = F, method = "UKF", eps = 1e-2,
                                  save_data = F, save_risk_set = F),
                   id = sims$res$id,
                   verbose = F,
                   max_T = 10)

  res <- do.call(ddhazard, arg_list)

  # matplot(sims$betas, type = "l")
  # matplot(res$state_vecs, add = T, type = "l")
  res <- res[c("state_vars", "state_vec", "Q")]
  # save_to_test(res, "UKF2")

  expect_equal(res, read_to_test("UKF2"))
})

test_that("Altering UKF alpha, beta and kappa change the results",{
  arg_list <- list(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                   by = 1,
                   data = sims$res,
                   Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                   Q = diag(rep(1e-1, ncol(sims$res) + 1 - 4)),
                   control = list(kappa = 01, alpha = 1, beta = 0,
                                  est_Q_0 = F, method = "UKF", eps = 1e-2),
                   id = sims$res$id,
                   verbose = F,
                   max_T = 10)
  suppressMessages(m1 <- do.call(ddhazard, arg_list))

  arg_list$control$beta <- 2

  suppressMessages(m2 <- do.call(ddhazard, arg_list))

  expect_true(class(all.equal(m1$state_vecs, m2$state_vecs)) == "character")
  expect_true(class(all.equal(m1$state_vars, m2$state_vars)) == "character")

  arg_list$control$alpha <- .2
  arg_list$control$kappa <- 4 / (.1)^2 + 4

  suppressMessages(m3 <- do.call(ddhazard, arg_list))

  expect_true(class(all.equal(m2$state_vecs, m3$state_vecs)) == "character")
  expect_true(class(all.equal(m2$state_vars, m3$state_vars)) == "character")

  arg_list$control$beta <- 0
  arg_list$control$alpha <- 1
  arg_list$control$kappa <- .5

  suppressMessages(m4 <- do.call(ddhazard, arg_list))

  expect_true(class(all.equal(m3$state_vecs, m4$state_vecs)) == "character")
  expect_true(class(all.equal(m3$state_vars, m4$state_vars)) == "character")
})

test_that("UKF on simulated data works with exponential model where both variables are used", {
  set.seed(9997)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 10, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = (1:10 - 5) / 2.5,
                            intercept_start = -5, sds = c(.1, rep(1, 10)))

  suppressMessages(result_exp <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(1, 11),
    Q = diag(1e-1, 11),
    control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3, method = "UKF",
                   debug = F, beta = 0, save_data = F, save_risk_set = F),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F,
    model = "exp_combined"))

  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, type = "l", lty = 2, add = T)

  result_exp <- result_exp[c("state_vars", "state_vec", "Q")]
  # save_to_test(result_exp, "UKF3")

  expect_equal(result_exp, read_to_test("UKF3"))
})


test_that("UKF on simulated data works with exponential models with only one of the variables", {
  set.seed(9997)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 10, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = (1:10 - 5) / 2.5,
                            intercept_start = -5, sds = c(.1, rep(1, 10)))

  suppressMessages(result_exp <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(1, 11),
    Q = diag(1e-1, 11),
    control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3, method = "UKF",
                   debug = F, beta = 0, save_data = F, save_risk_set = F),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F,
    model = "exp_bin"))

  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, type = "l", lty = 2, add = T)
  result_exp <- result_exp[c("state_vars", "state_vec", "Q")]
  # save_to_test(result_exp, "UKF4")

  expect_equal(result_exp, read_to_test("UKF4"))

  suppressMessages(result_exp <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(1, 11),
    Q = diag(1e-1, 11),
    control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3, method = "UKF",
                   debug = F, beta = 0, save_data = F, save_risk_set = F),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F,
    model = "exp_clip_time"))

  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, type = "l", lty = 2, add = T)
  result_exp <- result_exp[c("state_vars", "state_vec", "Q")]
  # save_to_test(result_exp, "UKF5")

  expect_equal(result_exp, read_to_test("UKF5"))

  suppressMessages(result_exp <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(1, 11),
    Q = diag(1e-1, 11),
    control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3, method = "UKF",
                   debug = F, beta = 0, save_data = F, save_risk_set = F),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F,
    model = "exp_clip_time_w_jump"))

  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, type = "l", lty = 2, add = T)
  result_exp <- result_exp[c("state_vars", "state_vec", "Q")]
  # save_to_test(result_exp, "UKF6")

  expect_equal(result_exp, read_to_test("UKF6"))
})

test_that("UKF second order model works", {
  for(m in c("logit", exp_model_names)){
    expect_no_error(result <- ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(1, 4),
      Q = diag(1e-1, 2),
      control = list(est_Q_0 = F, n_max = 10^4, eps = 10^-3,
                     method = "UKF", save_data = F, save_risk_set = F,
                     beta = 0),
      max_T = 30,
      id = head_neck_cancer$id, order = 2,
      verbose = F,
      model = m
    ))

    set.seed(9999)
    sim_f <- if(m == "logit") test_sim_func_logit else test_sim_func_exp
    sims <- sim_f(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 20,
                  lambda = 1/4,
                  x_range = 1, x_mean = 0.5, re_draw = T, beta_start = rep(0, 3),
                  intercept_start = -3, sds = c(.1, rep(.25, 3)),
                  tstart_sampl_func = function(...) max(0, runif(1, min = -4, max = 20 - 1)))

    expect_no_error(suppressWarnings(result_sim <-
                      ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                               by = 1,
                               data = sims$res,
                               Q_0 = diag(c(rep(.25, ncol(sims$res) + 1 - 4), rep(.25, ncol(sims$res) + 1 - 4))),
                               Q = diag(rep(.1, ncol(sims$res) + 1 - 4)),
                               control = list(est_Q_0 = F, method = "UKF", eps = 1e-1,
                                              ridge_eps = 1e-3),
                               id = sims$res$id,
                               verbose = F, model = m,
                               order = 2,
                               max_T = 20)))
  }
})

#############
# Test ukf for exponential model
test_that("UKF and exp models give previous computed results for head_neck",{
  suppressMessages(result_exp <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    Q_0 = diag(1, 2),
    Q = diag(1e-1, 2),
    control = list(est_Q_0 = F,
                   method = "UKF", beta = 0, alpha = 1, save_risk_set = F,
                   save_data = F),
    max_T = 30,
    id = head_neck_cancer$id, order = 1,
    verbose = F,
    model = "exp_clip_time_w_jump"))

  # plot(result_exp)
  result_exp <- result_exp[c("state_vecs", "state_vars", "Q")]
  # save_to_test(result_exp, "UKF_exp_head_neck")

  expect_equal(result_exp, read_to_test("UKF_exp_head_neck"))
})

# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
