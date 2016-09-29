# library(testthat); library(survival); library(parallel); source("R/test_utils.R")

result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2), Q_0 = diag(.1, 2),
  save_risk_set = T
)

for(use_parallel in c(T, F)){
  # Test that we get same estimate in one period estimates
  test_that(paste0("Testing one period in sample with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(1, 60)),
                       use_parallel = use_parallel)
    tmp_a_t_d_s <- rbind(result$a_t_d_s[-1, ], result$a_t_d_s[60, ])
    expect_equal(predict_$fits, c(result$hazard_func(tmp_a_t_d_s %*% c(1, 1))),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(0, 60)),
                       use_parallel = use_parallel)
    expect_equal(predict_$fits, c(result$hazard_func(tmp_a_t_d_s %*% c(1, 0))),
                 use.names = F, check.attributes = F)
  })

  # Check with two period estimates
  test_that(paste0("Testing two period in sample with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 2*(0:29), stop = 2*(1:30), group = rep(1, 30)),
                       use_parallel = use_parallel)
    tmp_a_t_d_s <- rbind(result$a_t_d_s[-1, ], result$a_t_d_s[60, ])
    fac1 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30) - 1, ] %*% c(1, 1)))
    fac2 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30), ] %*% c(1, 1)))
    expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 2*(0:29), stop = 2*(1:30), group = rep(0, 30)),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30) - 1, ] %*% c(1, 0)))
    fac2 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30), ] %*% c(1, 0)))
    expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
                 use.names = F, check.attributes = F)
  })

  # Check forcasting
  test_that(paste0("Testing forcasting with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 59, stop = 69, group = 1),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(result$a_t_d_s[60, ] %*% c(1, 1)))
    expect_equal(predict_$fits, (1 - (1 - fac1)^10),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 59, stop = 69, group = 0),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(result$a_t_d_s[60, ] %*% c(1, 0)))
    expect_equal(predict_$fits, (1 - (1 - fac1)^10),
                 use.names = F, check.attributes = F)
  })
}

# Check the terms
test_that("Testing term prediction",{
  for(g in 0:1){
    predict_ = predict(result, new_data = data.frame(start = 0:58, stop = 1:59, group = g))
    predict_terms = predict(result, new_data = data.frame(group = g), type = "term")
    expect_equal(dim(predict_terms$terms), c(60, 1, 2))
    predict_terms$terms <- predict_terms$terms[-1, , , drop = F]

    predict_terms_fit = exp(c(apply(predict_terms$terms, 1:2, sum))) /
      (exp(c(apply(predict_terms$terms, 1:2, sum))) + 1)
    expect_equal(unname(predict_$fits), predict_terms_fit, check.attributes = F)
  }
})

test_that("Terms prediction when parsing two observations",{
  predict_terms = predict(result, new_data = data.frame(group = 0:1), type = "term")
  expect_equal(dim(predict_terms$terms), c(60, 2, 2))
})

# Check for second order
result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2 * 2), Q_0 = diag(10, 2 * 2),
  Q = diag(c(1.0e-4, 1.0e-4, 0, 0)), est_Q_0 = F,
  n_max = 1e3,
  save_risk_set = T, order_ = 2
)

test_that("Calls with second order models do not throw errors", {
  for(g in c(0, 1))
    for(use_parallel in c(T, F))
      expect_no_error(bquote(
        predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(.(g), 60)),
                use_parallel = .(use_parallel))))
})

dum <- structure(list(model = "exponential"), "class" = class(result))

test_that("predict functions throws error when model is exponential",{
  expect_error(predict(dum))
})

warning("Implement test for second order random walk and predict")

# PBC dataset described in Fleming & Harrington (1991)
library(timereg)
data(pbc)
head(pbc)
# status at endpoint, 0/1/2 for censored, transplant, dead

pbc <- pbc[complete.cases(pbc[, c("time", "status", "age", "edema", "bili", "protime")]), ]
max(pbc$time[pbc$status == 2])
fit <- ddhazard(
  formula = survival::Surv(rep(0, nrow(pbc)), time, status == 2) ~ splines::ns(log(bili),df = 4),
  data = pbc, Q_0 = diag(rep(1e3, 5)), by = 100,
  Q = diag(rep(1e-2, 5)), max_T = 3600, est_Q_0 = F,
  verbose = F, save_risk_set = T)

predict(fit, new_data = pbc[1:5, ], type = "term")


######
# Exponential model
test_that("Test exponential model", expect_true(FALSE))
