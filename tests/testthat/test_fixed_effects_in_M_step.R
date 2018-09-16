context("Testing test_fixed_effects_in_M_step")

set.seed(548237)
sims <- test_sim_func_logit(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -4, sds = c(.1, rep(1, 3)))

test_that("Only fixed effects yields same results as bigglm with logit model", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    ddFixed_intercept() + ddFixed(x1) +
                    ddFixed(x2) + ddFixed(x3))

  suppressWarnings(
    res1 <- ddhazard(
      form, data = sims$res, model = "logit", by = 1, id = sims$res$id,
      max_T = 10, control = ddhazard_control(
        eps_fixed_parems = 1e-3, fixed_effect_chunk_size = 1e3,
        max_it_fixed_params = 10, fixed_terms_method = "M_step")))

  tmp_design <- get_survival_case_weights_and_data(
    form, data = sims$res, by = 1, id = sims$res$id,
    use_weights = F, max_T = 10)

  did_fail <- tryCatch(
    suppressWarnings(res2 <- bigglm(
      update(form, Y ~ . - ddFixed_intercept()), data = tmp_design$X,
      family = binomial(), chunksize = 1e3)),
    error = function(...) TRUE)
  if(isTRUE(did_fail))
    skip("bigglm failed (likely biglm error in 0.9-1 release)")

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects)))
})


test_that("Get previous results with logit model with some fixed terms", {
  dd_fit <- ddhazard(
    Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() +
      ddFixed(age) + ddFixed(log(albumin)) + edema +ddFixed(log(protime)) +
      log(bili),
    pbc2, id = pbc2$id, by = 100, max_T = 3600,
    Q_0 = diag(rep(100000, 2)), Q = diag(rep(0.001, 2)))

  # plot(dd_fit)
  dd_fit <- dd_fit[c("state_vars", "state_vecs","Q","fixed_effects")]
  # save_to_test(dd_fit, "fixed_terms_pbc")

  expect_equal(dd_fit, read_to_test("fixed_terms_pbc"))
})


test_that("Only fixed effects yields same results as bigglm and static_glm with exponential model", {
  sims <- exp_sim_200

  form <- formula(survival::Surv(tstart, tstop, event) ~
                    ddFixed_intercept() + ddFixed(x1) + ddFixed(x2) +
                    ddFixed(x3))

  suppressWarnings(res1 <- ddhazard(
    form, data = sims$res, model = "exponential", by = 1,
    id = sims$res$id, max_T = 10, control = ddhazard_control(
      eps_fixed_parems = 1e-4, fixed_effect_chunk_size = 1e3,
      fixed_terms_method = "M_step")))

  tmp_design <- get_survival_case_weights_and_data(
    form, data = sims$res, by = 1, id = sims$res$id,
    use_weights = F, max_T = 10, is_for_discrete_model = F)

  did_fail <- tryCatch(
    res2 <- bigglm(
      update(form, Y ~ . - ddFixed_intercept() +
               offset(log(pmin(tstop, 10) - tstart))),
      data = tmp_design$X, family = poisson(), chunksize = 1e3,
      tolerance = 1e-4),
    error = function(...) TRUE)
  if(isTRUE(did_fail))
    skip("bigglm failed (likely biglm error in 0.9-1 release)")

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects))
               , tolerance = 1e-05)

  res3 <- static_glm(
    form, data = sims$res, max_T = 10, id = sims$res$id,
    family = "exponential")$coefficients

  expect_equal(res3, coef(res2))
})

test_that("Changing fixed effect control parems changes the result", {
  sims <- exp_sim_200
  cl <- quote(ddhazard(
    formula(survival::Surv(tstart, tstop, event) ~
              ddFixed_intercept() + ddFixed(x1) + ddFixed(x2) + ddFixed(x3)),
    data = sims$res, model = "exp_clip_time_w_jump", by = 1, id = sims$res$id,
    max_T = 10, control = ddhazard_control(
      eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3,
      fixed_terms_method = "M_step")))

  suppressWarnings(res1 <- eval(cl))

  # Should not make a difference
  cl_tmp <- cl
  ctrl <- eval(cl$control)
  ctrl$fixed_effect_chunk_size <- 1e4
  cl_tmp$control <- ctrl
  suppressWarnings(res2 <- eval(cl_tmp))
  expect_equal(res2$fixed_effects, res1$fixed_effects)

  # Should make a difference
  ctrl <- eval(cl$control)
  ctrl$eps_fixed_parems <- 1
  cl_tmp$control <- ctrl
  suppressWarnings(res2 <- eval(cl_tmp))
  expect_true(!all(res2$fixed_effects == res1$fixed_effects))

  # Should make a difference
  ctrl <- eval(cl$control)
  ctrl$max_it_fixed_params <- 5
  cl_tmp$control <- ctrl
  suppressWarnings(res2 <- eval(cl_tmp))
  expect_true(!all(res2$fixed_effects == res1$fixed_effects))
})

test_that("Gets previous results with exponential model with some fixed terms", {
  sims <- exp_sim_500
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    ddFixed_intercept() + x1 + x2 + x3)

  res1 <- ddhazard(
    form, data = sims$res, model = "exponential", by = 1,
    id = sims$res$id, max_T = 10, control = ddhazard_control(
      eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3,
      save_risk_set = F, save_data = F, n_max = 1e2,
      fixed_terms_method = "M_step"),
    Q_0 = diag(100000, 3), Q = diag(.01, 3))

  # matplot(sims$betas[, 1:4], type = "l", lty = 2
  #         ylim = range(sims$betas, res1$state_vecs, res1$fixed_effects))
  # matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
  # abline(h = res1$fixed_effects)
  # res1$fixed_effects
  res1 <- res1[c("state_vars","state_vecs","fixed_effects")]
  expect_known_value(res1, "fixed_terms1.RDS", update = FALSE)
})

test_that("UKF with fixed effects works", {
  sims <- exp_sim_200

  fit <- ddhazard(
    formula(survival::Surv(tstart, tstop, event) ~
              ddFixed_intercept() + ddFixed(x1) + x2 + x3),
    Q = diag(1, 2), Q_0 = diag(10, 2),
    data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10,
    control = ddhazard_control(
      method = "UKF", fixed_parems_start = rep(0, 2),
      save_data = F, save_risk_set = F,
      fixed_terms_method = "M_step"))

  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(fit$state_vecs, type = "l", lty = 2, col = 3:4, add = T)
  # abline(h = fit$fixed_effects, col = 1:2, lty = 2)
  fit <- fit[c("state_vars", "state_vecs", "fixed_effects")]
  # save_to_test(fit, "fixed_terms_UKF")

  expect_equal(fit, read_to_test("fixed_terms_UKF"))

  #####
  fit <- ddhazard(formula(survival::Surv(tstart, tstop, event) ~
                            ddFixed_intercept() + ddFixed(x1) + x2 + x3),
                  Q = diag(.1, 2), Q_0 = diag(10, 2),
                  data = sims$res, model = "exponential",
                  by = 1, id = sims$res$id, max_T = 10,
                  control = ddhazard_control(
                    method = "UKF", fixed_terms_method = "M_step"))


  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(fit$state_vecs, type = "l", lty = 2, col = 3:4, add = T)
  # abline(h = fit$fixed_effects, col = 1:2, lty = 2)

  fit <- fit[c("state_vars", "state_vecs", "fixed_effects")]
  expect_known_value(fit, "fixed_terms_UKF_exp.RDS", update = FALSE)
})


test_that("posterior_approx gives previous found values with fixed effects in M-step", {
  skip_on_cran()

  set.seed(950466)
  f1 <- suppressWarnings(ddhazard(
    Surv(tstart, tstop, death == 2) ~ ddFixed(age) + ddFixed(edema) +
      log(albumin) + log(protime) + log(bili), pbc2,
     id = pbc2$id, by = 100, max_T = 3600,
     control = ddhazard_control(
       method = "SMA",  fixed_terms_method = "M_step",
       eps_fixed_parems = 1e-4),
     Q_0 = diag(rep(100000, 4)), Q = diag(rep(1e-3, 4))))

  # plot(f1)
  # f1$fixed_effects

  f1 <- f1[c("state_vecs", "state_vecs")]
  expect_known_value(f1, "posterior_approx_logit_fixed_M.RDS", update = FALSE)
})


test_that("Only fixed effects yields same results as bigglm and static_glm with exponential model with weights", {
  sims <- exp_sim_200
  set.seed(2555647)
  ws <- runif(nrow(sims$res))
  ws <- ws * (nrow(sims$res) / sum(ws))
  sims$res <- cbind(sims$res, ws = ws)

  form <- formula(survival::Surv(tstart, tstop, event) ~
                    ddFixed_intercept() + ddFixed(x1) + ddFixed(x2) +
                    ddFixed(x3))

  suppressWarnings(res1 <- ddhazard(
    form, data = sims$res, model = "exponential", by = 1,
    control = ddhazard_control(
      eps_fixed_parems = 1e-4, fixed_effect_chunk_size = 1e3,
      fixed_terms_method = "M_step"),
    weights = sims$res$ws,  id = sims$res$id, max_T = 10))

  tmp_design <- get_survival_case_weights_and_data(
    form, data = sims$res, by = 1, id = sims$res$id,
    use_weights = F, max_T = 10, is_for_discrete_model = F)

  did_fail <- tryCatch(
    res2 <- bigglm(
      update(form, Y ~ . -ddFixed_intercept() +
               offset(log(pmin(tstop, 10) - tstart))),
             data = tmp_design$X, family = poisson(), chunksize = 1e3,
             tolerance = 1e-4, weights = ~ ws),
    error = function(...) TRUE)
  if(isTRUE(did_fail))
    skip("bigglm failed (likely biglm error in 0.9-1 release)")

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects))
               , tolerance = 1e-05)

  res3 <- static_glm(
    form, data = sims$res, max_T = 10, id = sims$res$id,
    family = "exponential", weights = sims$res$ws)$coefficients

  expect_equal(res3, coef(res2))
})
