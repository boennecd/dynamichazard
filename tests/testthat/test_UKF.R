context("Testing UKF")

test_that("UKF throws error when first sigm points weight is zero", {
  expect_error({
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(1, 2), a_0 = c(-3, 0),
      Q = diag(1e-1, 2),
      control = ddhazard_control(kappa = 0, method = "UKF"),
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
    control = ddhazard_control(
      eps = 10^-3, method = "UKF", beta = 0, alpha = 1),
    max_T = 30,
    id = head_neck_cancer$id, order = 1,
    verbose = F,
    model = "logit"
  ))

  # plot(result)
  result <- result[c("state_vecs", "state_vars", "lag_one_cov")]
  # save_to_test(result, "UKF1")

  expect_equal(result, read_to_test("UKF1"))
})

test_that("UKF does not fail and both methods give the same",{
  sims <- logit_sim_200
  res_new <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id,
    by = 1,
    data = sims$res,
    a_0 = rep(0, ncol(sims$res) + 1 - 4),
    Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
    Q = diag(rep(1, ncol(sims$res) + 1 - 4)),
    verbose = F,
    control = ddhazard_control(eps = 1e-2, est_Q_0 = F, method = "UKF",
                               kappa = 0.1, alpha = 1, beta = 0),
    id = sims$res$id,
    max_T = 10)

  res_old <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id,
    by = 1,
    data = sims$res,
    a_0 = rep(0, ncol(sims$res) + 1 - 4),
    Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
    Q = diag(rep(1, ncol(sims$res) + 1 - 4)),
    control = ddhazard_control(
      est_Q_0 = F, kappa = 0.1, alpha = 1, beta = 0, eps = 1e-2, method = "UKF_org"),
    id = sims$res$id,
    max_T = 10)

  expect_equal(res_new$state_vecs, res_old$state_vecs)
  expect_equal(res_new$state_vars, res_old$state_vars)
})

test_that("Changing time scale in UKF does no change results when other parems are changed accoridngly",{
  sims <- logit_sim_200

  arg_list <- list(formula = survival::Surv(tstart, tstop, event) ~ . - id,
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
  sims <- logit_sim_200

  arg_list <- list(formula = survival::Surv(tstart, tstop, event) ~ . - id,
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

  # matplot(sims$betas, type = "l", lty = 1, ylim = range(sims$betas, res$state_vecs))
  # matplot(res$state_vecs, add = T, type = "l", lty = 2)
  res <- res[c("state_vars", "state_vec", "Q")]
  # save_to_test(res, "UKF2")

  expect_equal(res, read_to_test("UKF2"))
})

test_that("Altering UKF alpha, beta and kappa change the results",{
  sims <- logit_sim_200

  arg_list <- list(formula = survival::Surv(tstart, tstop, event) ~ . - id,
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

test_that("UKF works on simulated data works with exponential model and gives previous results", {
  sims <- exp_sim_200

  result_exp <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(1, 11),
    Q = diag(1e-1, 11),
    control = ddhazard_control(
      eps = 10^-2, method = "UKF",
      debug = F, beta = 0, save_data = F, save_risk_set = F),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F,
    model = "exponential")

  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, type = "l", lty = 2, add = T)
  result_exp <- result_exp[c("state_vars", "state_vecs", "Q")]
  expect_known_value(result_exp, "UKF4.RDS", update = FALSE)
})

test_that("UKF second order model works (that is, gives no errors...)", {
  # TODO: make a better test. Change the description if you do
  for(m in c("logit", exp_model_names)){
    expect_no_error(result <- ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(1, 4),
      Q = diag(1e-1, 2),
      control = ddhazard_control(est_Q_0 = F, method = "UKF", beta = 0),
      max_T = 30,
      id = head_neck_cancer$id, order = 2,
      verbose = F,
      model = m
    ))

    sims <- if(m == "logit") logit_sim_200 else exp_sim_200
    expect_no_error(
      result_sim <-
        ddhazard(
          formula = survival::Surv(tstart, tstop, event) ~ . - id,
           by = 1,
           data = sims$res,
           Q_0 = diag(c(
             rep(.1,  ncol(sims$res) + 1 - 4),
             rep(.1, ncol(sims$res) + 1 - 4))),
           Q = diag(.01, ncol(sims$res) + 1 - 4),
           control = ddhazard_control(
             est_Q_0 = F, method = "UKF",
             eps = .1 # Just want see the a few iterations passes
           ),
           id = sims$res$id,
           model = m,
           order = 2,
           max_T = 10))
  }
})
