# Test first order
test_that("Expecting plot calls to succed with first order model", {
  result = ddhazard(
    formula = Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    order_ = 1
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
    by = 1, control = list(eps = 1e-2, est_Q_0 = F),
    a_0 = rep(0, 4), Q_0 = diag(1, 4),
    Q = diag(c(5e-3, 5e-3, 0, 0)),
    order_ = 2,
    max_T = 40
  )

  expect_no_error(plot(result, type = "cov", cov_index = 1))
  expect_no_error(plot(result, type = "cov", cov_index = 2))
})

warning("Implement more plot tests for other than EKF logit model")
