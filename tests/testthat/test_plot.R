context("Testing plot functions")

# Test first order
test_that("Expecting plot calls to succed with first order model", {
  cl <- quote(ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(1e-3, 2),
    order = 1))

  for(i in 0:1){
    if(i) cl$control <- quote(ddhazard_control(method = "UKF"))

    result = eval(cl)

    if(i < 1){
      space_error <- residuals(result, type ="std_space_error")
      expect_no_error(plot(space_error, result))
    }

    expect_no_error(plot(result, type = "cov", cov_index = 1))
    expect_no_error(plot(result, type = "cov", cov_index = 2))
    expect_no_error(plot(result, type = "cov"))
  }

  suppressMessages(pbc_fit <- ddhazard(
    formula = survival::Surv(tstart/100, tstop/100, death == 2) ~ log(bili) + log(protime),
    data = pbc2, model = "exp_clip_time", by = 1, max_T = 36,
    Q_0 = diag(2, 3), Q = diag(1e-3, 3), verbose = F,
    id = pbc2$id,
    control = ddhazard_control(LR = 1, eps = 1e-3, save_risk_set = F)))

  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 1))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 2))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 3))
  expect_no_error(plot(pbc_fit, type = "cov"))
})

# Test second order
test_that("Expecting plot calls to succed with second order model", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, control = ddhazard_control(eps = 1e-2, est_Q_0 = F),
    a_0 = rep(0, 4), Q_0 = diag(1, 4),
    Q = diag(c(5e-3, 5e-3)),
    order = 2,
    max_T = 40)

  expect_no_error(plot(result, type = "cov", cov_index = 1))
  expect_no_error(plot(result, type = "cov", cov_index = 2))
  expect_no_error(plot(result, type = "cov"))

  suppressMessages(pbc_fit <- ddhazard(
    formula = survival::Surv(tstart/100, tstop/100, death == 2) ~ log(bili) + log(protime),
    data = pbc2, model = "exp_clip_time", by = 1, max_T = 36,
    Q_0 = diag(5, 6), Q = diag(c(rep(1e-3, 3))),
    id = pbc2$id, order = 2,
    control = ddhazard_control(LR = .01, eps = 1e-2, save_risk_set = F)))

  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 1))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 2))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 3))
  expect_no_error(plot(pbc_fit, type = "cov"))
})

test_that("Alters mfcol and sets it back", {
  set.seed(747)
  sims <- test_sim_func_exp(n_series = 2e2, n_vars = 10, t_0 = 0, t_max = 10,
                    x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                    intercept_start = -5, sds = c(.1, rep(1, 10)))

  suppressMessages(result_exp <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(10, 11),
    Q = diag(1e-2, 11),
    control = ddhazard_control(
      est_Q_0 = F, eps = 10^-2, n_max = 10^3,
      save_data = F, save_risk_set = F, denom_term = 1e-2),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F,
    model = "exp_clip_time"))

  for(i in 1:10){
    expect_no_error(plot(result_exp, type = "cov", cov_index = 1:i))
    expect_equal(getOption("mfcol"), NULL)
  }
})
