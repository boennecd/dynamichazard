# Simple test that methods calls succeds
result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2), Q_0 = diag(.1, 2),
  save_risk_set = T
)

test_that("Calls to residuals should succed",{
  expect_no_error(residuals(result, "stdSpaceErrors"))
  expect_no_error(residuals(result, type = "Pearson", data_ = head_neck_cancer))
  expect_no_error(residuals(result, type = "Raw", data_ = head_neck_cancer))
})

dum <- structure(list(model = "poisson"), "class" = class(result))

test_that("residuals functions throws error when model is poisson",{
  expect_error(predict(dum))
})

warning("Implement more test for residuals method (e.g. second order models, the actual values etc.)")
