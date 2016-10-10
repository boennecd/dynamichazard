# library(testthat); library(survival); library(parallel); source("R/test_utils.R")

result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2), Q_0 = diag(.1, 2))

for(use_parallel in c(T, F)){
  # Test that we get same estimate in one period estimates
  test_that(paste0("Testing one period in sample with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(1, 60)),
                       use_parallel = use_parallel)
    tmp_state_vecs <- rbind(result$state_vecs[-1, ], result$state_vecs[60, ])
    expect_equal(predict_$fits, c(result$hazard_func(tmp_state_vecs %*% c(1, 1))),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(0, 60)),
                       use_parallel = use_parallel)
    expect_equal(predict_$fits, c(result$hazard_func(tmp_state_vecs %*% c(1, 0))),
                 use.names = F, check.attributes = F)
  })

  # Check with two period estimates
  test_that(paste0("Testing two period in sample with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 2*(0:29), stop = 2*(1:30), group = rep(1, 30)),
                       use_parallel = use_parallel)
    tmp_state_vecs <- rbind(result$state_vecs[-1, ], result$state_vecs[60, ])
    fac1 = c(result$hazard_func(tmp_state_vecs[2*(1:30) - 1, ] %*% c(1, 1)))
    fac2 = c(result$hazard_func(tmp_state_vecs[2*(1:30), ] %*% c(1, 1)))
    expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 2*(0:29), stop = 2*(1:30), group = rep(0, 30)),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(tmp_state_vecs[2*(1:30) - 1, ] %*% c(1, 0)))
    fac2 = c(result$hazard_func(tmp_state_vecs[2*(1:30), ] %*% c(1, 0)))
    expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
                 use.names = F, check.attributes = F)
  })

  # Check forcasting
  test_that(paste0("Testing forcasting with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 59, stop = 69, group = 1),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(result$state_vecs[60, ] %*% c(1, 1)))
    expect_equal(predict_$fits, (1 - (1 - fac1)^10),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 59, stop = 69, group = 0),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(result$state_vecs[60, ] %*% c(1, 0)))
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
  Q = diag(c(1.0e-4, 1.0e-4, 0, 0)),
  control = list(n_max = 1e3, save_risk_set = T, est_Q_0 = F),
  order_ = 2
)

test_that("Calls with second order models do not throw errors", {
  for(g in c(0, 1))
    for(use_parallel in c(T, F))
      expect_no_error(bquote(
        predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(.(g), 60)),
                use_parallel = .(use_parallel))))
})

dum <- structure(list(model = "exponential"), "class" = class(result))

test_that("predict functions throws error when model it is exponential",{
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
  Q = diag(rep(1e-2, 5)), max_T = 3600,
  control = list(est_Q_0 = F))

predict(fit, new_data = pbc[1:5, ], type = "term")


######
# Exponential model
fit <- ddhazard(
  formula = survival::Surv(tstart, tstop, status == 2) ~
    age + log(bili) + log(protime),
  data = pbc2, Q_0 = diag(rep(1e3, 4)), by = 100,
  Q = diag(rep(1e-2, 4)), max_T = 3600,
  model = "exponential", control = list(est_Q_0 = F, LR = .4))


test_that("Terms from predict with exponential outcome are correct", {
  pred <- predict(fit, new_data = pbc2, type = "term", tstart = "tstart", tstop = "tstop")

  expect_equal(dim(pred$terms), c(1 + 3600/100, nrow(pbc2), 4))
  tmp_mat <- t(cbind(rep(1, nrow(pbc2)), pbc2$age, log(pbc2$bili), log(pbc2$protime)))

  for(i in 1:ncol(fit$state_vecs))
    expect_equal(pred$terms[, , i], fit$state_vecs[, i] %o% tmp_mat[i,])

  expect_message(
    respone_pred <- predict(fit, new_data = pbc2, type = "response", tstart = "tstart", tstop = "tstop"),
    "start and stop times \\('tstart' and 'tstop'\\) are in data. Prediction will match these periods")

  set.seed(192301258)
  rand_indicies <- sample.int(nrow(pbc2), 1000)
  test_rows <- pbc2[rand_indicies, ]
  test_rows <- cbind(rep(1, nrow(test_rows)),
                     test_rows$age, log(test_rows$bili), log(test_rows$protime),
                     test_rows)

  for(j in seq_along(test_rows)){
    tmp <- test_rows[j, ]
    bins_breaks <- seq(0, 3600, by = 100)
    start <- findInterval(tmp$tstart, bins_breaks)
    stop  <- findInterval(tmp$tstop, bins_breaks, left.open = T)
    p_survival <- 1

    t_min <- max(bins_breaks[start], tmp$tstart)
    for(i in start:stop){
      t_max <- min(bins_breaks[i + 1], tmp$tstop)
      p_survival <- p_survival * exp(- exp(fit$state_vecs[i + 1, ] %*% unlist(tmp[, 1:4])) * (t_max - t_min))
      t_min <- t_max
    }

    expect_equal(unname(respone_pred$fits[rand_indicies[j]]), c(1 - p_survival), info = "")
  }
})
