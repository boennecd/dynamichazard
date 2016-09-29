###############
# Simple test that methods calls succeds
result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1, max_T = 40,
  a_0 = rep(0, 2), Q_0 = diag(10, 2),
  Q = diag(1e-2, 2),
  save_risk_set = T
)

test_that("Calls to residuals should succed",{
  expect_no_error(residuals(result, "std_space_error"))
  expect_no_error(residuals(result, "space_error"))
  expect_no_error(residuals(result, type = "pearson", data_ = head_neck_cancer))
  expect_no_error(residuals(result, type = "raw", data_ = head_neck_cancer))
})

dum <- structure(list(model = "exponential"), "class" = class(result))
test_that("residuals functions throws error when model is exponential",{
  expect_error(predict(dum))
})


########
# Test state space errors

test_that("State space error gives previous result with logit model", {
  std_res <- residuals(result, "std_space_error")

  expect_true(std_res$standardize)
  expect_s3_class(std_res, "fahrmeier_94_SpaceErrors")

  non_std <- residuals(result, "space_error")
  expect_true(!non_std$standardize)
  expect_s3_class(non_std, "fahrmeier_94_SpaceErrors")

  expect_equal(non_std$Covariances, std_res$Covariances)

  for(i in seq_len(nrow(std_res$residuals)))
    expect_equal(unname(std_res$residuals[i, ]),
                 c(solve(t(chol(non_std$Covariances[, , i]))) %*% non_std$residuals[i, ]))

})

# PBC dataset described in Fleming & Harrington (1991)
library(timereg)
data(pbc)
head(pbc)
# status at endpoint, 0/1/2 for censored, transplant, dead

# pbc <- pbc[complete.cases(pbc[, c("time", "status", "age", "edema", "bili", "protime")]), ]
# max(pbc$time[pbc$status == 2])
# fit <- ddhazard(
#   formula = survival::Surv(rep(0, nrow(pbc)), time, status == 2) ~
#     age + edema + log(bili) + log(protime),
#   data = pbc, Q_0 = diag(rep(1e3, 5)), by = 100,
#   Q = diag(rep(1e-2, 5)), max_T = 3600, est_Q_0 = F,
#   verbose = F)
#
# resids <- residuals(fit, "std_space_error")
# matplot(resids$residuals)
# resids <- residuals(fit, "space_error")
# matplot(resids$residuals)

warning("Implement more test for residuals method (e.g. second order models, the actual values etc.)")
