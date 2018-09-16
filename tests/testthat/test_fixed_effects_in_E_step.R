context("Testing fixed_effects_in_E_step")

test_that("Works with EKF and logit model and with one fixed and one non-fixed", {
          fit <- ddhazard(
            survival::Surv(start, stop, event) ~ ddFixed(group),
            Q_0 = matrix(1), Q = matrix(.1),
            data = head_neck_cancer,
            by = 1,
            control = ddhazard_control(
              fixed_terms_method = "E_step", save_risk_set = F, save_data = F))

          # plot(fit)
          # fit$fixed_effects
          fit <- fit[c("state_vecs","state_vars","Q","fixed_effects")]
          # save_to_test(fit, "E_step_fixed_head_neck")

          expect_equal(fit, read_to_test("E_step_fixed_head_neck"), tolerance = 1.490116e-08)
})

test_that("Gives previous results with call as for pbc as in vignette", {
  dd_fit <- ddhazard(
    Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() +
      ddFixed(age) + ddFixed(log(albumin)) + edema +ddFixed(log(protime)) + log(bili),
    pbc2, id = pbc2$id, by = 100, max_T = 3600,
    Q_0 = diag(rep(100000, 2)), Q = diag(rep(0.001, 2)),
    control = ddhazard_control(fixed_terms_method = "E_step"))

  # plot(dd_fit)
  # dd_fit$fixed_effects

  dd_fit <- dd_fit[c("state_vars", "state_vecs", "Q", "fixed_effects")]
  # save_to_test(dd_fit, "E_step_fixed_pbc")

  expect_equal(dd_fit, read_to_test("E_step_fixed_pbc"), tolerance = 1.490116e-08)
})

test_that("Works with EKF and continous time model and that it work with one fixed and one non-fixed", {
  fit <- ddhazard(survival::Surv(start, stop, event) ~ ddFixed(group),
                  Q_0 = matrix(100000), Q = matrix(.1),
                  data = head_neck_cancer, max_T = 35,
                  by = 1, model = "exponential",
                  control = ddhazard_control(
                    fixed_terms_method = "E_step", save_risk_set = F,
                    save_data = F))

  # plot(fit$state_vecs)
  # fit$fixed_effects
  fit <- fit[c("state_vecs", "state_vars","Q","fixed_effects")]
  expect_known_value(fit, "E_step_fixed_exp_head_neck.RDS", tolerance = 1e-8,
                     update = FALSE)
})

set.seed(849239)
sims_logit <- test_sim_func_logit(n_series = 1e3, n_vars = 4, x_range = 2, t_max = 10, x_mean = 0,
                                  beta_start = rnorm(4), intercept_start = -4)
# sum(sims_logit$res$event)
form <- formula(survival::Surv(tstart, tstop, event) ~ ddFixed_intercept() + ddFixed(x1) + x2 + x3 + x4)

test_that("Works with UKF and logit model", {
 fit <- ddhazard(form, Q_0 = diag(1, 3), Q = diag(.1, 3),
                  data = sims_logit$res, id = sims_logit$res$id,
                  by = 1,
                  control = ddhazard_control(
                    fixed_terms_method = "E_step",
                    save_risk_set = F, save_data = F,
                    method = "UKF"))


  # matplot(sims_logit$betas, lty = 1, type = "l")
  # matplot(fit$state_vecs, lty = 2, type = "l", add = T, col = 3:5)
  # abline(h = fit$fixed_effects, col = 1:2, lty = 2)

  fit <- fit[c("state_vars", "state_vecs", "Q", "lag_one_cov")]
  # save_to_test(fit, "E_step_sim_UKF")

  expect_equal(fit, read_to_test("E_step_sim_UKF"), tolerance = 1.490116e-08)
})

test_that("Works with second order random walk and logit model", {
  fit <- ddhazard(form, Q_0 = diag(c(rep(1000000, 3), rep(10, 3))), Q = diag(.01, 3),
                  data = sims_logit$res, id = sims_logit$res$id,
                  by = 1, order = 2,
                  control = ddhazard_control(
                    fixed_terms_method = "E_step",
                    save_risk_set = F, save_data = F, LR = .5))

  # matplot(sims_logit$betas, lty = 1, type = "l")
  # matplot(fit$state_vecs[, 1:3], lty = 2, type = "l", add = T, col = 3:5)
  # abline(h = fit$fixed_effects, col = 1:2, lty = 2)

  fit <- fit[c("state_vars", "state_vecs", "Q", "lag_one_cov")]
  # save_to_test(fit, "E_step_sim_EKF_order_two")

  expect_equal(fit, read_to_test("E_step_sim_EKF_order_two"), tolerance = 0.1)
})

test_that("Get loglike to work so you can use verbose with E-step fixed effects", {
  # There was a bug - thus this test is added
  expect_output(
    fit <- ddhazard(form, Q_0 = diag(c(rep(1, 3), rep(.1, 3))), Q = diag(1e-4, 3),
                    data = sims_logit$res, id = sims_logit$res$id,
                    by = 1, order = 2,
                    control = ddhazard_control(
                      fixed_terms_method = "E_step", save_risk_set = F,
                      save_data = F, LR = 1 / 1.25),
                    verbose = 5)
    , regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood of the mean path is\\s+")
})



set.seed(1534834)
sims_exp <- test_sim_func_exp(n_series = 4e2, n_vars = 3, x_range = 1, t_max = 10, x_mean = 0,
                              beta_start = 1, intercept_start = -3, is_fixed = 1:2)
# sum(sims_exp$res$event)

form <- formula(survival::Surv(tstart, tstop, event) ~ ddFixed_intercept() + ddFixed(x1) + x2 + x3)

test_that("Works with UKF and continous time model", {
  fit <- ddhazard(form, Q_0 = diag(1, 2), Q = diag(1, 2),
                  data = sims_exp$res, id = sims_exp$res$id,
                  by = 1, model = "exponential", max_T = 10,
                  control = ddhazard_control(
                    fixed_terms_method = "E_step", save_risk_set = F,
                    save_data = F, method = "UKF"))

  # matplot(sims_exp$betas, lty = 1, type = "l",
  #         ylim = range(fit$fixed_effects, fit$state_vecs, sims_exp$betas))
  # matplot(fit$state_vecs, lty = 2, type = "l", add = T, col = 3:4)
  # abline(h = fit$fixed_effects, col = 1:2, lty = 2)

  fit <- fit[c("state_vars", "state_vecs", "Q", "lag_one_cov")]
  expect_known_value(fit, "E_step_sim_exp_UKF.RDS", update = FALSE)
})

test_that("Works with second order random walk and continous time model",{
  fit <- ddhazard(form,
                  Q_0 = diag(1, 4), Q = diag(1, 2),
                  data = sims_exp$res, id = sims_exp$res$id,
                  by = 1, model = "exp_clip_time_w_jump", max_T = 10,
                  control = ddhazard_control(
                    fixed_terms_method = "E_step", save_risk_set = F,
                    save_data = F, method = "UKF"),
                  order = 2)

  # matplot(sims_exp$betas, lty = 1, type = "l",
  #         ylim = range(fit$fixed_effects, fit$state_vecs, sims_exp$betas))
  # matplot(fit$state_vecs[, 1:2], lty = 2, type = "l", add = T, col = 3:4)
  # abline(h = fit$fixed_effects, col = 1:2, lty = 2)

  fit <- fit[c("state_vars", "state_vecs", "Q", "lag_one_cov")]
  expect_known_value(fit, "E_step_sim_exp_UKF_order_two.RDS", update = FALSE)
})

test_that("posterior_approx gives previous found values with fixed effects in E-step", {
  set.seed(950466)
  f1 <- ddhazard(Surv(tstart, tstop, death == 2) ~ ddFixed(age) + ddFixed(edema) +
                   log(albumin) + log(protime) + log(bili), pbc2,
                 id = pbc2$id, by = 100, max_T = 3600,
                 control = ddhazard_control(
                   method = "SMA",  fixed_terms_method = "E_step"),
                 Q_0 = diag(rep(100000, 4)), Q = diag(rep(0.01, 4)))

  # plot(f1)
  # f1$fixed_effects

  f1 <- f1[c("state_vecs", "state_vecs")]
  expect_known_value(f1, "posterior_approx_logit_fixed_E.RDS", update = FALSE)
})
