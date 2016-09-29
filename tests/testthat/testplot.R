# Test first order
result = ddhazard(
  formula = Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2), Q_0 = diag(1, 2),
  save_risk_set = T, order_ = 1
)

test_that("Update plot tests", expect_true(FALSE))

# space_error <- residuals(result, type ="stdSpaceErrors")
#
# test_that("Expecting the following calls to succed with a second order fit", {
#   expect_no_error(plot(space_error, result))
#
#   expect_no_error(plot(result, type = "cov", cov_index = 1))
#   expect_no_error(plot(result, type = "cov", cov_index = 2))
#
#   expect_no_error(plot(result, new_data = data.frame(group = 1), method = "delta"))
#   expect_no_error(plot(result, new_data = data.frame(group = 0), method = "delta"))
#   expect_no_error(plot(result, new_data = data.frame(group = 1), method = "backtransform"))
#   expect_no_error(plot(result, new_data = data.frame(group = 0), method = "backtransform"))
# })
#
# # Test second order
# result = ddhazard(
#   formula = Surv(start, stop, event) ~ group,
#   data = head_neck_cancer,
#   by = 1,
#   a_0 = rep(0, 2 * 2), Q_0 = diag(10, 2 * 2),
#   Q = diag(c(1.0e-4, 1.0e-4, 0, 0)), est_Q_0 = F,
#   n_max = 1e3,
#   save_risk_set = T, order_ = 2
# )
#
# test_that("Expecting the following calls to succed with a second order fit", {
#   expect_no_error(plot(result, type = "cov", cov_index = 1))
#   expect_no_error(plot(result, type = "cov", cov_index = 2))
#
#   expect_no_error(plot(result, new_data = data.frame(group = 1), method = "delta"))
#   expect_no_error(plot(result, new_data = data.frame(group = 0), method = "delta"))
#   expect_no_error(plot(result, new_data = data.frame(group = 1), method = "backtransform"))
#   expect_no_error(plot(result, new_data = data.frame(group = 0), method = "backtransform"))
# })
#
# dum <- structure(list(model = "poisson"), "class" = class(result))
#
# test_that("Plot functions throws error when model is poisson",{
#   expect_error(plot(dum))
# })

warning("Implement more plot tests")
