if(interactive()){
  library(survival); library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")
}


tmp <- file("tmp.txt")
sink(tmp)
result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  control = list(
    est_Q_0 = F, save_data = F, save_risk_set = F,
    method = "posterior_approx", debug = T),
  Q_0 = diag(1e0, 2), Q = diag(0.01, 2),
  max_T = 45,
  id = head_neck_cancer$id, order = 1)
sink()
close(tmp)

plot(result)

test_that("Taylor function for logit works",
          expect_true(FALSE))
