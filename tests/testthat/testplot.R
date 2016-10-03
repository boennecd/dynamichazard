# Test first order
test_that("Expecting plot calls to succed with first order model", {
  arg_list <- list(
    formula = Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    order_ = 1)

  for(i in 0:1){
    if(i) arg_list$control <- list("method" = "UKF")

    result = do.call(ddhazard, arg_list)

    if(i < 1){
      space_error <- residuals(result, type ="std_space_error")
      expect_no_error(plot(space_error, result))
    }

    expect_no_error(plot(result, type = "cov", cov_index = 1))
    expect_no_error(plot(result, type = "cov", cov_index = 2))
  }

  pbc_fit <- ddhazard(
    formula = survival::Surv(tstart, tstop, status == 2) ~ log(bili) + log(protime),
    data = pbc2, model = "exponential", by = 100, max_T = 3600,
    Q_0 = diag(10, 3), Q = diag(1e-3, 3), verbose = F,
    id = pbc2$id,
    control = list(LR = .5, eps = 1e-3, save_risk_set = F))

  space_error <- residuals(pbc_fit, type ="std_space_error")
  expect_no_error(plot(space_error, pbc_fit))

  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 1))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 2))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 3))
})

# Test second order
test_that("Expecting plot calls to succed with second order model", {
  result = ddhazard(
    formula = Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, control = list(eps = 1e-2, est_Q_0 = F),
    a_0 = rep(0, 4), Q_0 = diag(1, 4),
    Q = diag(c(5e-3, 5e-3, 0, 0)),
    order_ = 2,
    max_T = 40
  )

  expect_no_error(plot(result, type = "cov", cov_index = 1))
  expect_no_error(plot(result, type = "cov", cov_index = 2))

  pbc_fit <- ddhazard(
    formula = survival::Surv(tstart, tstop, status == 2) ~ log(bili) + log(protime),
    data = pbc2, model = "exponential", by = 100, max_T = 3600,
    Q_0 = diag(10, 6), Q = diag(c(rep(1e-3, 3), rep(0, 3))),
    id = pbc2$id, order_ = 2,
    control = list(LR = .35, eps = 1e-2, save_risk_set = F))

  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 1))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 2))
  expect_no_error(plot(pbc_fit, type = "cov", cov_index = 3))
})

warning("Implement more plot tests for other UKF and exponential fits")
