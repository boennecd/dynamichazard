# Test first order
test_that("Expecting plot calls to succed with first order model", {
  result = ddhazard(
    formula = Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    save_risk_set = T, order_ = 1
  )

  space_error <- residuals(result, type ="std_space_error")
  expect_no_error(plot(space_error, result))

  expect_no_error(plot(result, type = "cov", cov_index = 1))
  expect_no_error(plot(result, type = "cov", cov_index = 2))
})

# Test second order
test_that("Expecting plot calls to succed with first order model", {
  result = ddhazard(
    formula = Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, eps = 1e-2,
    a_0 = rep(0, 4), Q_0 = diag(1, 4),
    Q = diag(c(5e-3, 5e-3, 0, 0)),
    save_risk_set = T, order_ = 2,
    est_Q_0 = F, max_T = 40
  )

  expect_no_error(plot(result, type = "cov", cov_index = 1))
  expect_no_error(plot(result, type = "cov", cov_index = 2))
})

warning("Implement more plot tests for other than EKF logit model")
