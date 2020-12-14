context("Testing predict")

result <- ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2),
  Q_0 = diag(.1, 2),
  Q = diag(1, 2))

test_that("Testing one period in sample", {
  predict_ = predict(result, new_data = data.frame(
    start = 0:59, stop = 1:60, group = factor(rep(1, 60))))
  tmp_state_vecs <- rbind(result$state_vecs[-1, ], result$state_vecs[60, ])
  expect_equal(
    predict_$fits, c(result$discrete_hazard_func(
      tmp_state_vecs %*% c(1, 1), 1)),
    use.names = F, check.attributes = F)

  predict_ = predict(
    result, new_data = data.frame(start = 0:59, stop = 1:60,
                                  group = factor(rep(2, 60))))
  expect_equal(
    predict_$fits, c(result$discrete_hazard_func(tmp_state_vecs %*% c(1, 0), 1)),
    use.names = F, check.attributes = F)
})

test_that("Testing two period in sample", {
  predict_ = predict(
    result, new_data = data.frame(start = 2*(0:29), stop = 2*(1:30),
                                  group = factor(rep(1, 30))))
  tmp_state_vecs <- rbind(result$state_vecs[-1, ], result$state_vecs[60, ])
  fac1 = c(result$discrete_hazard_func(
    tmp_state_vecs[2*(1:30) - 1, ] %*% c(1, 1), 1))
  fac2 = c(result$discrete_hazard_func(
    tmp_state_vecs[2*(1:30), ] %*% c(1, 1), 1))
  expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
               use.names = F, check.attributes = F)

  predict_ = predict(result, new_data = data.frame(
    start = 2*(0:29), stop = 2*(1:30), group = factor(rep(2, 30))))
  fac1 = c(result$discrete_hazard_func(
    tmp_state_vecs[2*(1:30) - 1, ] %*% c(1, 0), 1))
  fac2 = c(result$discrete_hazard_func(
    tmp_state_vecs[2*(1:30), ] %*% c(1, 0), 1))
  expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
               use.names = F, check.attributes = F)
})

test_that("Testing forcasting with predict", {
  predict_ = predict(
    result, new_data = data.frame(start = 59, stop = 69, group = factor(1)))
  fac1 = c(result$discrete_hazard_func(
    result$state_vecs[60, ] %*% c(1, 1), 1))
  expect_equal(predict_$fits, (1 - (1 - fac1)^10),
               use.names = F, check.attributes = F)

  predict_ = predict(result, new_data = data.frame(
    start = 59, stop = 69, group = factor(2)))
  fac1 = c(result$discrete_hazard_func(
    result$state_vecs[60, ] %*% c(1, 0), 1))
  expect_equal(predict_$fits, (1 - (1 - fac1)^10),
               use.names = F, check.attributes = F)
})

# Check the terms
test_that("Testing term prediction",{
  for(g in 1:2){
    predict_ = predict(
      result, new_data = data.frame(start = 0:58, stop = 1:59,
                                    group = factor(g)))
    predict_terms = predict(result, new_data = data.frame(group = factor(g)),
                            type = "term")
    expect_equal(dim(predict_terms$terms[[1L]]), c(59, 2))
    predict_terms$terms <- predict_terms$terms[[1L]]
    predict_terms_fit <- 1 / (1 + exp(-rowSums(predict_terms$terms)))
    expect_equal(unname(predict_$fits), predict_terms_fit, check.attributes = F)
  }
})

test_that("Term prediction with fixed effects", {
  fit <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1, a_0 = 0, Q_0 = as.matrix(1), Q = 1)

  for(g in 2:1){
    predict_terms = predict(
      fit, new_data = data.frame(group = as.factor(g)), type = "term")

    expect_equal(c(predict_terms$terms[[1L]]), c(unname(fit$state_vecs))[-1L],
                 check.attributes = F)
    expect_true(all(predict_terms$fixed_terms[[1L]] ==
                      (g == 1) * unname(fit$fixed_effects)))
  }

  suppressWarnings(fit <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ -1 + ddFixed(group),
    data = head_neck_cancer, by = 1))

  for(g in 1:2){
    predict_terms = predict(
      fit, new_data = data.frame(group = factor(x = g)),
      type = "term")

    expect_equal(dim(predict_terms$terms[[1L]]), c(59, 0))
    .which <- grepl(paste0(g, "$"), names(fit$fixed_effects))
    expect_true(all(
      predict_terms$fixed_terms[[1L]] == unname(fit$fixed_effects[.which])))
  }
})

test_that("Terms prediction when parsing two observations",{
  predict_terms = predict(result, new_data = data.frame(group = factor(2:1)),
                          type = "term")
  expect_length(predict_terms, 6L)
  expect_length(predict_terms$terms, 2)
  expect_length(predict_terms$fixed_terms, 2)
  expect_length(predict_terms$tstart, 2)
  expect_length(predict_terms$tstop, 2)

  expect_equal(dim(predict_terms$terms[[1L]]), c(59, 2))
  expect_equal(dim(predict_terms$terms[[2L]]), c(59, 2))

  expect_length(predict_terms$fixed_terms[[1L]], 59)
  expect_length(predict_terms$fixed_terms[[2L]], 59)

  expect_length(predict_terms$tstart[[1L]], 59)
  expect_length(predict_terms$tstart[[2L]], 59)

  expect_length(predict_terms$tstop[[1L]], 59)
  expect_length(predict_terms$tstop[[2L]], 59)
})

test_that("Testing one period in sample with fixed effects", {
  .data <- head_neck_cancer
  .data$group <- relevel(.data$group, ref = 2)

  result <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = .data,
    by = 1, a_0 = 0, Q_0 = as.matrix(1), Q = 1)

  for(g in 1:2){
    predict_ = predict(result, new_data =
                         data.frame(start = 0:59, stop = 1:60,
                                    group = factor(rep(g, 60))))
    tmp_state_vecs <-
      rbind(result$state_vecs[-1, , drop = F], result$state_vecs[60, ])
    expect_equal(
      as.vector(predict_$fits),
      c(result$discrete_hazard_func(
        tmp_state_vecs + result$fixed_effects * (g - 1), 1)),
      use.names = F, check.attributes = F)
  }

  # ddhazard throws a warning when static_glm could be used instead
  result <- suppressWarnings(ddhazard(
    formula = survival::Surv(start, stop, event) ~ -1 + ddFixed(group),
    data = .data, by = 1))

  for(g in 1:2){
    predict_ = predict(
      result, new_data = data.frame(
        start = 0:59, stop = 1:60, group = factor(rep(g, 60))))
    tmp_state_vecs <- rbind(
      result$state_vecs[-1, , drop = F], result$state_vecs[60, ])
    expect_equal(
      unname(predict_$fits),
      rep(result$discrete_hazard_func(
        c(c(g == 1, g == 2) %*% result$fixed_effects), 1), 60),
      use.names = F, check.attributes = F)
  }
})


test_that("a warning is provided with predict if max_threads or use_parallel is used", {
  # check that we get a warning but the result is the same
  wa <- "'max_threads' and 'use_parallel' is not used anymore"
  for(ty in c("term", "response")){
    p1 <- predict(result, new_data = head_neck_cancer[1, ], type = "term")
    expect_warning(
      p2 <- predict(result, new_data = head_neck_cancer[1, ], type = "term",
                    use_parallel = TRUE),
      regexp = wa, fixed = TRUE)
    expect_warning(
      p3 <- predict(result, new_data = head_neck_cancer[1, ], type = "term",
                    max_threads = TRUE),
      regexp = wa, fixed = TRUE)
    expect_warning(
      p4 <- predict(result, new_data = head_neck_cancer[1, ], type = "term",
                    use_parallel = TRUE, max_threads = TRUE),
      regexp = wa, fixed = TRUE)
    expect_equal(p1, p2)
    expect_equal(p1, p3)
    expect_equal(p1, p4)
  }
})

test_that("predict yields the correct covariance matrix", {
  fit <- suppressWarnings(ddhazard(
    Surv(tstart, tstop, event) ~ cbind(x1, x1^2, x1^3),
    data = logit_sim_200$res, max_T = 10, model = "logit",
    control = ddhazard_control(n_max = 1), Q = diag(1, 4),
    Q_0 = diag(1, 4), by = 1))

  t_data <- within(
    logit_sim_200$res[1:2, ], {
      tstart <- floor(tstart)
      tstop <- ceiling(tstop)
    })
  preds <- predict(
    fit, new_data = t_data, tstart = "tstart", tstop = "tstop",
    type = "term", sds = TRUE)

  for(i in 1:nrow(t_data)){
    td <- t_data[i, ]
    XJ <- matrix(
      c(1,     0,       0,       0,
        0, td$x1, td$x1^2, td$x1^3), ncol = 2)
    pd <- preds$varcov[[i]]
    sa <- findInterval(td$tstart, fit$times, left.open = FALSE)
    so <- findInterval(td$tstop , fit$times, left.open = TRUE)
    i <- 0L
    for(j in sa:so){
      i <- i + 1L
      expect_equal(
        pd[, , i], crossprod(XJ, fit$state_vars[, , j + 1L] %*% XJ))
    }
  }

  # just run this to make sure that is does not fail when we ask for future
  # sds
  t_data$tstop[2] <- fit$times[length(fit$times)] + 5L
  preds <- predict(
    fit, new_data = t_data, tstart = "tstart", tstop = "tstop",
    type = "term", sds = TRUE)

  tis <- fit$times[1L]:(fit$times[length(fit$times)] + 5L)
  cv_mats <- array(dim = c(dim(fit$state_vars)[1:2], length(tis)))
  cv_mats[, , seq_along(fit$times)] <- fit$state_vars
  for(i in (length(fit$times) + 1L):dim(cv_mats)[3])
    cv_mats[, , i] <- cv_mats[, , i - 1L] + fit$Q

  for(i in 1:nrow(t_data)){
    td <- t_data[i, ]
    XJ <- matrix(
      c(1,     0,       0,       0,
        0, td$x1, td$x1^2, td$x1^3), ncol = 2)
    pd <- preds$varcov[[i]]
    sa <- findInterval(td$tstart, tis, left.open = FALSE)
    so <- findInterval(td$tstop , tis, left.open = TRUE)
    i <- 0L
    for(j in sa:so){
      i <- i + 1L
      expect_equal(
        pd[, , i], crossprod(XJ, cv_mats[, , j + 1L] %*% XJ))
    }
  }

  # has correct dimension when only one variable is time-varying
  fit <- suppressWarnings(ddhazard(
    Surv(tstart, tstop, event) ~ ddFixed(cbind(x1, x1^2, x1^3)),
    data = logit_sim_200$res, max_T = 10, model = "logit",
    control = ddhazard_control(n_max = 1), Q = 1, Q_0 = 1, by = 1))
  preds <- predict(
    fit, new_data = t_data[1, ], tstart = "tstart", tstop = "tstop",
    type = "term", sds = TRUE)
  expect_equal(dim(preds$sds[[1]]), dim(preds$terms[[1L]]))
  expect_equal(dim(preds$varcov[[1]])[3:2], dim(preds$terms[[1L]]))
  expect_equal(dim(preds$varcov[[1]])[1], dim(preds$terms[[1L]])[2])
})



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
                    Q = diag(1, 6), Q_0 = diag(10, 6),
                    control = ddhazard_control(LR = .5)))

  s$res$tstart_ceil <- ceiling(s$res$tstart)
  s$res$tstop_ceil <- as.integer(pmin(ceiling(s$res$tstop), 10))

  suppressMessages(preds <- predict(fit, new_data = s$res, tstart = "tstart_ceil", tstop = "tstop_ceil"))
  preds <- 1 - tapply(preds$fits, s$res$id, function(x) prod(1 - x))

  other_s <- get_survival_case_weights_and_data(
    formula = survival::Surv(tstart, tstop, event) ~ . - id,
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
  # tmp_id <- 10
  # dumdum[dumdum[, "id"] == tmp_id, ]
  # s$res[s$res$id == tmp_id, ]
  # other_s[other_s$id == tmp_id, ]
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
  # Check for second order
  result <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 2 * 2), Q_0 = diag(10, 2 * 2),
    Q = diag(c(1.0e-4, 1.0e-4)),
    control = ddhazard_control(n_max = 1e3, save_risk_set = T, est_Q_0 = F),
    order = 2)

  for(g in c(0, 1))
    for(use_parallel in c(T, F))
      expect_no_error(bquote(
        predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(.(g), 60)),
                use_parallel = .(use_parallel))))
})

test_that("predict functions throws error when one of the old exponential models are used",{
  dum <- structure(list(model = "exp_combined"), "class" = class(result))
  expect_error(predict(dum), regexp = "Functions for model 'exp_combined' is not implemented",
               fixed = TRUE)
})

test_that("Term predict work with second order random walk", {
  f1 <- ddhazard(
    formula = Surv(tstart, tstop, death == 2) ~ log(bili),
    data = pbc2,
    order = 2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = ddhazard_control(method = "GMA", LR = .33),
    Q_0 = diag(1, 4), Q = diag(rep(1e-4, 2)))

  dum <- data.frame(
    bili = exp(1))

  preds <- predict(f1, new_data = dum, type = "term", sds = T)

  # plot(f1)
  expect_equal(preds$terms[[1L]], f1$state_vecs[-1, 1:2],
               check.attributes = FALSE)
  expect_equal(preds$sds[[1L]],
               sqrt(t(apply(f1$state_vars, 3, diag))[-1, 1:2]),
               check.attributes = FALSE)
})

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
  model = "exponential"))

# plot(fit)

test_that("Terms from predict with exponential outcome are correct", {
  skip_on_cran()

  pred <- predict(
    fit, new_data = pbc2_l, type = "term", tstart = "tstart", tstop = "tstop")

  expect_length(pred$terms, nrow(pbc2_l))
  tmp_mat <- with(pbc2_l, t(cbind(1, age, log(bili), log(protime))))
  for(j in 1:nrow(pbc2_l)){
    trs <- pred$terms[[j]]
    idx <- findInterval(pred$tstart[[j]], fit$times,
                        rightmost.closed = FALSE) + 1L
    idx <- idx[idx <= nrow(fit$state_vecs)]
    keep <- seq_along(idx)

    expect_equal(
      t(trs[keep, ]), apply(
        fit$state_vecs[idx, , drop = FALSE], 1, "*", y = tmp_mat[, j]),
      check.attributes = FALSE)
  }

  expect_message(
    respone_pred <- predict(
      fit, new_data = pbc2_l, type = "response",
      tstart = "tstart", tstop = "tstop"),
    "start and stop times \\('tstart' and 'tstop'\\) are in data. Prediction will match these periods")

  set.seed(192301258)
  rand_indicies <- sample.int(nrow(pbc2_l), 100)
  test_rows <- pbc2_l[rand_indicies, ]
  test_rows <- cbind(rep(1, nrow(test_rows)),
                     test_rows$age, log(test_rows$bili), log(test_rows$protime),
                     test_rows)

  for(j in seq_along(test_rows)){
    tmp <- test_rows[j, ]
    bins_breaks <- fit$times
    start <- findInterval(tmp$tstart, bins_breaks)
    stop  <- findInterval(tmp$tstop , bins_breaks, left.open = T)
    p_survival <- 1

    t_min <- max(bins_breaks[start], tmp$tstart)
    for(i in start:stop){
      t_max <- min(bins_breaks[i + 1], tmp$tstop)
      if(t_max == t_min)
        break
      fa <- exp(
        - exp(fit$state_vecs[i + 1, ] %*%
                unlist(tmp[, 1:4])) * (t_max - t_min))
      p_survival <- p_survival * fa
      t_min <- t_max
    }

    expect_equal(unname(respone_pred$fits[rand_indicies[j]]),
                 c(1 - p_survival))
  }
})




test_that("Different non-integer time_scales gives the same result with predict results", {
  skip_on_cran()
  scales <- exp(seq(-2, 2, .05))

  fit_exp <- expression(suppressMessages(local({
    fit <- ddhazard(
      formula = survival::Surv(start * .by, stop * .by, event) ~ group,
      data = head_neck_cancer,
      by = .by,
      control = ddhazard_control(est_Q_0 = F, save_data = F),
      a_0 = rep(0, 2), Q_0 = diag(1e2, 2), Q = diag(1e-2 / .by, 2),
      id = head_neck_cancer$id, order = 1)

    tmp_dat <- head_neck_cancer[1, ]
    tmp_dat$stop <- 40 * .by
    p1 <- predict(fit, new_data = tmp_dat)
    t1 <- predict(fit, new_data = tmp_dat, type = "term")

    tmp_dat <- head_neck_cancer[1, ]
    tmp_dat$start <- 20 * .by
    tmp_dat$stop <- 100 * .by
    p2 <- predict(fit, new_data = tmp_dat)
    t2 <- predict(fit, new_data = tmp_dat, type = "term")

    list(p1 = p1$fits, p2 = p2$fits, t1 = t1$terms, t2 = t2$terms)
    })))

  .by <- scales[1]
  f1 <- eval(fit_exp)
  for(.by in scales[-1]){
    f2 <- eval(fit_exp)
    expect_equal(f1, f2, info = paste0("by = ", .by))
  }
})
