context("Testing testfixed_terms_est_in_M_step")

set.seed(548237)
sims <- test_sim_func_logit(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -4, sds = c(.1, rep(1, 3)))

test_that("Only fixed effects yields same results as bigglm with logit model", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  suppressWarnings(
    res1 <- ddhazard(form, data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10,
                     control = list(eps_fixed_parems = 1e-3, fixed_effect_chunk_size = 1e3, max_it_fixed_params = 10,
                                    fixed_terms_method = "M_step")))

  tmp_design <- get_survival_case_weights_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10)

  suppressWarnings(res2 <- bigglm(update(form, Y ~ .), data = tmp_design$X, family = binomial(), chunksize = 1e3))

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects)))
})


test_that("Get previous results with logit model with some fixed terms", {
  dd_fit <- ddhazard(
    Surv(tstart, tstop, death == 2) ~ ddFixed(1) +
      ddFixed(age) + ddFixed(log(albumin)) + edema +ddFixed(log(protime)) + log(bili),
    pbc2, id = pbc2$id, by = 100, max_T = 3600,
    Q_0 = diag(rep(100000, 2)), Q = diag(rep(0.001, 2)))

  # plot(dd_fit)
  dd_fit <- dd_fit[c("state_vars", "state_vecs","Q","fixed_effects")]
  # save_to_test(dd_fit, "fixed_terms_pbc")

  expect_equal(dd_fit, read_to_test("fixed_terms_pbc"))
})


test_that("Only fixed effects yields same results as bigglm with exponential model", {
  set.seed(4682146)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -4, sds = c(.1, rep(1, 3)))

  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  suppressWarnings(res1 <- ddhazard(form, data = sims$res, model = "exp_clip_time_w_jump", by = 1, id = sims$res$id, max_T = 10,
                                    control = list(eps_fixed_parems = 1e-4, fixed_effect_chunk_size = 1e3,
                                                   fixed_terms_method = "M_step")))

  tmp_design <- get_survival_case_weights_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10, is_for_discrete_model = F)

  suppressWarnings(res2 <- bigglm(update(form, Y ~ . + offset(log(pmin(tstop, t) - pmax(tstart, t - 1)))),
                                  data = tmp_design$X, family = poisson(), chunksize = 1e3,
                                  tolerance = 1e-4))

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects))
               , tolerance = 1e-05)
})

test_that("Changing fixed effect control parems changes the result", {
  arg_list <- list(
    formula(survival::Surv(tstart, tstop, event) ~
              -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3)),
    data = sims$res, model = "exp_clip_time_w_jump", by = 1, id = sims$res$id, max_T = 10,
    control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3,
                   fixed_terms_method = "M_step"))

  suppressWarnings(res1 <- do.call(ddhazard, arg_list))

  # Should not make a difference
  arg_list_tmp <- arg_list
  arg_list_tmp$control$fixed_effect_chunk_size <- 1e4
  suppressWarnings(res2 <- do.call(ddhazard, arg_list_tmp))
  expect_equal(res2$fixed_effects, res1$fixed_effects)

  # Should make a difference
  arg_list_tmp <- arg_list
  arg_list_tmp$control$eps_fixed_parems <- 1
  suppressWarnings(res2 <- do.call(ddhazard, arg_list_tmp))
  expect_true(!all(res2$fixed_effects == res1$fixed_effects))

  # Should make a difference
  arg_list_tmp <- arg_list
  arg_list_tmp$control$max_it_fixed_params <- 5
  suppressWarnings(res2 <- do.call(ddhazard, arg_list_tmp))
  expect_true(!all(res2$fixed_effects == res1$fixed_effects))
})

test_that("Get previous results with exponential model with some fixed terms", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + x1 + x2 + x3)

  suppressMessages(
    res1 <- ddhazard(form, data = sims$res, model = "exp_clip_time_w_jump", by = 1, id = sims$res$id, max_T = 10,
                    control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3,
                                   save_risk_set = F, save_data = F, n_max = 1e2,
                                   denom_term = .0001,
                                   fixed_terms_method = "M_step"),
                    Q_0 = diag(rep(100000, 3)), Q = diag(.01, 3)))

  # matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
  # matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
  # res1$fixed_effects
  res1 <- res1[c("state_vars","state_vecs","fixed_effects")]
  # save_to_test(res1, "fixed_terms1")

  expect_equal(res1, read_to_test("fixed_terms1"))
})

test_that("UKF with fixed effects works", {
  set.seed(2231412)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -4, sds = c(.1, rep(1, 3)))

  fit <- ddhazard(formula(survival::Surv(tstart, tstop, event) ~
                            -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + x2 + x3),
                  data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10,
                  control = list(method = "UKF", fixed_parems_start = rep(0, 2),
                                 save_data = F, save_risk_set = F,
                                 fixed_terms_method = "M_step"))


  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(fit$state_vecs, type = "l", lty = 2, col = 3:4, add = T)
  fit <- fit[c("state_vars", "state_vecs", "fixed_effects")]
  # save_to_test(fit, "fixed_terms_UKF")

  expect_equal(fit, read_to_test("fixed_terms_UKF"))

  #####
  fit <- ddhazard(formula(survival::Surv(tstart, tstop, event) ~
                            -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + x2 + x3),
                  data = sims$res, model = "exp_clip_time_w_jump", by = 1, id = sims$res$id, max_T = 10,
                  control = list(method = "UKF",
                                 fixed_terms_method = "M_step"))


  # matplot(sims$betas, type = "l", ylim = range(fit$state_vecs, sims$betas))
  # matplot(fit$state_vecs, type = "l", col = 3:4, add = T, lty = 1)

  fit <- fit[c("state_vars", "state_vecs", "fixed_effects")]
  # save_to_test(fit, "fixed_terms_UKF_exp")

  expect_equal(fit, read_to_test("fixed_terms_UKF_exp"))
})


test_that("posterior_approx gives previous found values with fixed effects in M-step", {
  set.seed(950466)
  f1 <- ddhazard(Surv(tstart, tstop, death == 2) ~ ddFixed(age) + ddFixed(edema) +
                  log(albumin) + log(protime) + log(bili), pbc2,
                 id = pbc2$id, by = 100, max_T = 3600,
                 control = list(method = "SMA",  fixed_terms_method = "M_step"),
                 Q_0 = diag(rep(100000, 4)), Q = diag(rep(0.01, 4)))

  # plot(f1)
  # f1$fixed_effects
  f1 <- f1[c("state_vecs", "state_vecs")]
  # save_to_test(f1, "posterior_approx_logit_fixed_M")

  expect_equal(f1, read_to_test("posterior_approx_logit_fixed_M"), tolerance = 1.490116e-08)
})


test_that("Only fixed effects yields same results as bigglm with exponential model with weights", {
  set.seed(2555647)
  ws <- runif(nrow(sims$res))
  ws <- ws * (nrow(sims$res) / sum(ws))
  sims$res <- cbind(sims$res, ws = ws)

  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  suppressWarnings(res1 <- ddhazard(form, data = sims$res, model = "exp_clip_time_w_jump", by = 1, id = sims$res$id, max_T = 10,
                                    control = list(eps_fixed_parems = 1e-4, fixed_effect_chunk_size = 1e3,
                                                   fixed_terms_method = "M_step"),
                                    weights = sims$res$ws))

  tmp_design <- get_survival_case_weights_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10, is_for_discrete_model = F)

  suppressWarnings(res2 <- bigglm(
    update(form, Y ~ . + offset(log(pmin(tstop, t) - pmax(tstart, t - 1)))),
           data = tmp_design$X, family = poisson(), chunksize = 1e3,
           tolerance = 1e-4, weights = ~ ws))

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects))
               , tolerance = 1e-05)
})
