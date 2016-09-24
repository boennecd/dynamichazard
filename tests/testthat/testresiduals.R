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

# PBC dataset described in Fleming & Harrington (1991)
library(timereg)
data(pbc)
head(pbc)
# status at endpoint, 0/1/2 for censored, transplant, dead


pbc <- pbc[complete.cases(pbc[, c("time", "status", "age", "edema", "bili", "protime")]), ]
max(pbc$time[pbc$status == 2])
fit <- ddhazard(
  formula = survival::Surv(rep(0, nrow(pbc)), time, status == 2) ~
    age + edema + log(bili) + log(protime),
  data = pbc, Q_0 = diag(rep(1e3, 5)), by = 100,
  Q = diag(rep(1e-2, 5)), max_T = 3600, est_Q_0 = F,
  verbose = T)

plot(fit, cov_index = 5, type = "cov")

warning("Implement more test for residuals method (e.g. second order models, the actual values etc.)")
