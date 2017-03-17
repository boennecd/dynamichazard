if(interactive()){
  rm(list = ls())
  library(testthat); library(survival); library(dynamichazard); library(parallel)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")
}

# Had issues with win builder. Thus, these lines
test_name <- "predict"
cat("\nRunning", test_name, "\n")

suppressMessages(result <- ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2), Q_0 = diag(.1, 2)))

for(use_parallel in c(T, F)){
  # Test that we get same estimate in one period estimates
  test_that(paste0("Testing one period in sample with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(1, 60)),
                       use_parallel = use_parallel, max_threads = 1)
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

test_that("Term prediction with fixed effects",{
  fit <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1, a_0 = 0, Q_0 = as.matrix(1))

  for(g in 0:1){
    predict_terms = predict(fit, new_data = data.frame(group = g), type = "term")

    expect_equal(c(predict_terms$terms[,,]), c(unname(fit$state_vecs)),
                 check.attributes = F)
    expect_equal(c(predict_terms$fixed_terms), c((g == 1) * unname(fit$fixed_effects)))
  }

  suppressWarnings(fit <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ -1 + ddFixed(group),
    data = head_neck_cancer, by = 1))

  for(g in 1:2){
    predict_terms = predict(fit, new_data = data.frame(group = factor(x = g, levels = 1:2)), type = "term")

    expect_equal(dim(predict_terms$terms), c(60, 1, 0))
    expect_equal(c(predict_terms$fixed_terms), unname(fit$fixed_effects[g]))
  }
})

test_that("Terms prediction when parsing two observations",{
  predict_terms = predict(result, new_data = data.frame(group = 0:1), type = "term")
  expect_equal(dim(predict_terms$terms), c(60, 2, 2))
})


for(use_parallel in c(T, F)){
  # Test that we get same estimate in one period estimates
  test_that(paste0("Testing one period in sample with fixed effects and use_parallel ", use_parallel),{
    result <- ddhazard(
      formula = survival::Surv(start, stop, event) ~ ddFixed(group),
      data = head_neck_cancer,
      by = 1, a_0 = 0, Q_0 = as.matrix(1))

    for(g in c(0, 1)){
      predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(g, 60)),
                         use_parallel = use_parallel)
      tmp_state_vecs <- rbind(result$state_vecs[-1, , drop = F], result$state_vecs[60, ])
      expect_equal(predict_$fits, c(result$hazard_func(tmp_state_vecs %*% 1 + c(result$fixed_effects * g))),
                   use.names = F, check.attributes = F)
    }

    suppressWarnings(result <- ddhazard(
      formula = survival::Surv(start, stop, event) ~ -1 + ddFixed(group),
      data = head_neck_cancer, by = 1))

    for(g in c(1, 2)){
      predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = factor(rep(g, 60), levels = 1:2)),
                         use_parallel = use_parallel)
      tmp_state_vecs <- rbind(result$state_vecs[-1, , drop = F], result$state_vecs[60, ])
      expect_equal(unname(predict_$fits), rep(result$hazard_func(c(c(g == 1, g == 2) %*% result$fixed_effects)), 60),
                   use.names = F, check.attributes = F)
    }
  })
}

# Check for second order
result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2 * 2), Q_0 = diag(10, 2 * 2),
  Q = diag(c(1.0e-4, 1.0e-4)),
  control = list(n_max = 1e3, save_risk_set = T, est_Q_0 = F),
  order = 2
)

test_that("get_survival_case_weights_and_data and predict yields consistent result with ids", {
  set.seed(1992)

  s <- test_sim_func_logit(
    n_series = 1e3,
    n_vars = 5,
    beta_start = c(-1, -.5, 0, 1.5, 2),
    intercept_start = -4,
    sds = c(.1, rep(1, 5)),
    t_max = 10,
    x_range = 1,
    x_mean = .5)

  suppressMessages(
    fit <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3 + x4 + x5,
                    data = s$res, max_T = 10, by = 1, id = s$res$id,
                    control = list(LR = .5)))

  s$res$tstart_ceil <- ceiling(s$res$tstart)
  s$res$tstop_ceil <- as.integer(pmin(ceiling(s$res$tstop), 10))

  suppressMessages(preds <- predict(fit, new_data = s$res, tstart = "tstart_ceil", tstop = "tstop_ceil"))
  preds <- 1 - tapply(preds$fits, s$res$id, function(x) prod(1 - x))

  other_s <- get_survival_case_weights_and_data(
    formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - id - event,
    data = s$res, max_T = 10, by = 1, id = s$res$id,
    use_weights = F)$X

  other_s$tstart <- other_s$t - 1
  other_s$tstop <- other_s$t

  suppressMessages(other_preds <- predict(fit, new_data = other_s, tstart = "tstart", tstop = "tstop"))
  other_preds <- 1 - tapply(other_preds$fits, other_s$id, function(x) prod(1 - x))

  expect_equal(other_preds, preds)

  # # May be usefull for debugging
  # unique(other_s$id)[abs(other_preds - preds) > 1e-14]
  # dumdum <- cbind(other_preds - preds, id = unique(other_s$id))[abs(other_preds - preds) > 1e-14,]
  #
  # dumdum[dumdum[, "id"] == 5, ]
  # s$res[s$res$id == 5, ]
  # other_s[other_s$id == 5, ]
  #
  # suppressMessages(preds <- predict(fit, new_data = s$res, tstart = "tstart_ceil", tstop = "tstop_ceil"))
  # preds$fits[s$res$id == 5]
  #
  # suppressMessages(other_preds <- predict(fit, new_data = other_s, tstart = "tstart", tstop = "tstop"))
  # other_preds$fits[other_s$id == 5]
  #
  # (1 - prod(1 - other_preds$fits[other_s$id == 5][1:2])) - preds$fits[s$res$id == 5][1]
  # 1 - prod(1 - other_preds$fits[other_s$id == 5][3:3]) - preds$fits[s$res$id == 5][2]
  # 1 - prod(1 - other_preds$fits[other_s$id == 5][4:5]) - preds$fits[s$res$id == 5][3]
  # 1 - prod(1 - other_preds$fits[other_s$id == 5][6:7]) - preds$fits[s$res$id == 5][4]
  # 1 - prod(1 - other_preds$fits[other_s$id == 5][8:9]) - preds$fits[s$res$id == 5][5]
  # 1 - prod(1 - other_preds$fits[other_s$id == 5][10:10]) - preds$fits[s$res$id == 5][6]
  #
  # d1 <- fit$state_vecs[9,] %*% c(1, unlist(s$res[s$res$id == 5, c("x1", "x2", "x3", "x4", "x5")][5, ]))
  # d2 <- fit$state_vecs[10,] %*% c(1, unlist(s$res[s$res$id == 5, c("x1", "x2", "x3", "x4", "x5")][5, ]))
  #
  # d1 <- exp(d1) / (1 + exp(d1))
  # d2 <- exp(d2) / (1 + exp(d2))
  #
  # d1 <- fit$state_vecs[10,] %*% c(1, unlist(s$res[s$res$id == 5, c("x1", "x2", "x3", "x4", "x5")][5, ]))
  # d2 <- fit$state_vecs[11,] %*% c(1, unlist(s$res[s$res$id == 5, c("x1", "x2", "x3", "x4", "x5")][5, ]))
  #
  # d1 <- exp(d1) / (1 + exp(d1))
  # d2 <- exp(d2) / (1 + exp(d2))
  #
  # 1 - (1 - d1) * (1 - d2)
})

test_that("Calls with second order models do not throw errors", {
  for(g in c(0, 1))
    for(use_parallel in c(T, F))
      expect_no_error(bquote(
        predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(.(g), 60)),
                use_parallel = .(use_parallel))))
})

dum <- structure(list(model = "exp_combined"), "class" = class(result))

test_that("predict functions throws error when model it is exponential",{
  expect_error(predict(dum))
})

# PBC dataset described in Fleming & Harrington (1991)
library(timereg, quietly = T)
data(pbc)
head(pbc)
# status at endpoint, 0/1/2 for censored, transplant, dead

pbc <- pbc[complete.cases(pbc[, c("time", "status", "age", "edema", "bili", "protime")]), ]
max(pbc$time[pbc$status == 2])
suppressMessages(fit <- ddhazard(
  formula = survival::Surv(rep(0, length(status)), time, status == 2) ~ splines::ns(log(bili),df = 4),
  data = pbc, Q_0 = diag(rep(1e3, 5)), by = 100,
  Q = diag(rep(1e-2, 5)), max_T = 3600,
  control = list(est_Q_0 = F, eps = .01)))

predict(fit, new_data = pbc[1:5, ], type = "term")

######
# Exponential model
pbc2_l <- pbc2
pbc2_l$tstart <- pbc2_l$tstart / 100
pbc2_l$tstop <- pbc2_l$tstop / 100

suppressMessages(fit <- ddhazard(
  formula = survival::Surv(tstart, tstop, death == 2) ~
    age + log(bili) + log(protime),
  data = pbc2_l, Q_0 = diag(rep(1e5, 4)), by = 1,
  id = pbc2_l$id,
  Q = diag(rep(1e-1, 4)), max_T = 36,
  model = "exp_clip_time_w_jump"))

# plot(fit)

test_that("Terms from predict with exponential outcome are correct", {
  pred <- predict(fit, new_data = pbc2_l, type = "term", tstart = "tstart", tstop = "tstop")

  expect_equal(dim(pred$terms), c(1 + 3600/100, nrow(pbc2_l), 4))
  tmp_mat <- t(cbind(rep(1, nrow(pbc2_l)), pbc2_l$age, log(pbc2_l$bili), log(pbc2_l$protime)))

  for(i in 1:ncol(fit$state_vecs))
    expect_equal(pred$terms[, , i], fit$state_vecs[, i] %o% tmp_mat[i,])

  expect_message(
    respone_pred <- predict(fit, new_data = pbc2_l, type = "response", tstart = "tstart", tstop = "tstop"),
    "start and stop times \\('tstart' and 'tstop'\\) are in data. Prediction will match these periods")

  set.seed(192301258)
  rand_indicies <- sample.int(nrow(pbc2_l), 100)
  test_rows <- pbc2_l[rand_indicies, ]
  test_rows <- cbind(rep(1, nrow(test_rows)),
                     test_rows$age, log(test_rows$bili), log(test_rows$protime),
                     test_rows)

  for(j in seq_along(test_rows)){
    tmp <- test_rows[j, ]
    bins_breaks <- seq(0, 36, by = 1)
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


# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
