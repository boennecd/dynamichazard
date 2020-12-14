context("Testing loglike")

test_that("ddhazard with verbose > 0 prints log likelihood",{
  for(method in c("EKF", "UKF")){
    for(o in c(1, 2)){
      Q_0_arg <- if(o == 1) diag(1, 2) else diag(1, 4)
      for(m in c("exponential", "logit")){
        if(m == "exponential" && method == "UKF")
          next
        eval(bquote(expect_output({
          ddhazard(survival::Surv(start, stop, event) ~ group,
                   data = head_neck_cancer, id = head_neck_cancer$id,
                   by = 1, max_T = 30, Q = diag(.01, 2),
                   Q_0 = .(Q_0_arg), model = .(m), order = .(o),
                   verbose = 5,
                   control = ddhazard_control(eps = .1, method = .(method)))
        }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood of the mean path is\\s+")))
}}}})

test_that("ddhazard with verbose > 0 prints log likelihood with fixed effects",{
  for(method in c("EKF", "UKF")){
    for(fixed_method in c("M_step", "E_step")){
      Q_0_arg <- if(method == "EKF") 1e5 else 1
      for(m in c("exponential", "logit")){
        if(m == "exponential" && method == "UKF")
          next
        eval(bquote(expect_output({
          ddhazard(survival::Surv(start, stop, event) ~ ddFixed(group),
                   data = head_neck_cancer, id = head_neck_cancer$id,
                   by = 2, max_T = 30, Q = diag(.1, 1),
                   Q_0 = .(Q_0_arg), model = .(m),
                   verbose = 5,
                   control = ddhazard_control(
                     eps = .1, method = .(method),
                     fixed_terms_method = .(fixed_method)))
        }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood of the mean path is\\s+")))
      }}}})

test_that("logLik for head_neck_cancer data set match previous results", {
  result <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(est_Q_0 = F),
    a_0 = rep(0, 2), Q_0 = diag(1e5, 2),
    Q = diag(1e-2, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  log_like <- logLik(object = result)
  result$data <- NULL
  result$risk_set <- NULL
  logLik(object = result, data = head_neck_cancer, id = head_neck_cancer$id)

  # plot(result)
  # print(log_like, digits = 16)

  old <- structure(-182.0072469079185,
                   class = "logLik")

  expect_equal(log_like, old)
})

test_that("Saving or not saving risk set or data gives the same result", {
  cl <- quote(ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    Q = diag(.2, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1,
    verbose = F))

  control_fit <- eval(cl)
  control_fit$logLik <- logLik(control_fit)

  for(save_risk_set in c(T, F))
    for(save_data in c(T, F)){
      if(save_risk_set && save_data)
        next

      cl$control <- ddhazard_control(
        save_risk_set = save_risk_set, save_data = save_data)
      new_fit <- eval(cl)
      new_fit$logLik <- logLik(
        new_fit, data = if(save_data) NULL else head_neck_cancer,
        id = if(save_risk_set) NULL else head_neck_cancer$id)

      expect_equal(new_fit$logLik, control_fit$logLik)
    }
})

test_that("logLik for head_neck_cancer data set with second order model", {
  result <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    Q_0 = diag(10, 4),
    Q = diag(1e-3, 2),
    control = ddhazard_control(est_Q_0 = F, eps = 2e-3),
    max_T = 45,
    id = head_neck_cancer$id, order = 2)

  log_like <- logLik(object = result, data_ = head_neck_cancer, id =  head_neck_cancer$id)


  # plot(result)
  # print(log_like, digits = 16)

  old <- structure(-59.32922017545229,
                   class = "logLik")

  expect_equal(log_like, old, tolerance = 1e-4)
})

test_that("logLik for head_neck_cancer data set match previous results with fixed effects", {
  result <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(
      est_Q_0 = F, fixed_terms_method = "M_step"),
    a_0 = 0, Q_0 = 1e5, Q = 1e-2,
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  # plot(result)
  log_like <- logLik(result, data = head_neck_cancer)

  # print(log_like, digits = 16)

  old <- structure(-254.7652723920431, class = "logLik")

  expect_equal(log_like, old, tolerance = 5e-3)
})

test_that("logLik for head_neck_cancer data with only fixed match bigglm", {
  form <- survival::Surv(start, stop, event) ~
    ddFixed_intercept() + ddFixed(group)

  suppressWarnings(result <- ddhazard(
    formula = form,
    data = head_neck_cancer,
    by = 1,
    max_T = 45,
    id = head_neck_cancer$id, order = 1,
    control = ddhazard_control(fixed_terms_method = "M_step")))

  tmp_design <- get_survival_case_weights_and_data(
    formula = form, data = head_neck_cancer, by = 1, max_T = 45, id = head_neck_cancer$id,
    use_weights = F)

  glm_fit <- glm(Y ~ as.factor(group), binomial(), tmp_design$X)

  tmp <- logLik(glm_fit)
  attributes(tmp)
  attr(tmp, "nobs") <- NULL
  attr(tmp, "df") <- NULL

  expect_equal(c(unname(result$fixed_effects)), unname(glm_fit$coefficients))
  expect_equal(logLik(result), tmp)
})


##############

test_that("logLik for simulated data versus old results", {
  sims <- logit_sim_200
  result = ddhazard(
    survival::Surv(tstart, tstop, event) ~ . - id,
    sims$res,
    by = 1, Q = diag(1e-2, 11),
    Q_0 = diag(1e5, 11),
    max_T = 10,
    id = sims$res$id, order = 1)

  # matplot(sims$betas, lty = 1, type = "l")
  # matplot(result$state_vecs, type = "l", lty = 2, add = T)

  log_like <- logLik(object = result)

  # dput(log_like)
  old <- structure(-136.706189501914,
                   class = "logLik")
  expect_equal(log_like, old)

  sims <- exp_sim_200
  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ . - id,
    sims$res,
    by = 1, Q = diag(1e-2, 11),
    Q_0 = diag(1e5, 11),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F, model = "exponential",
    control = ddhazard_control(n_max = 150))

  # matplot(sims$betas, lty = 1, type = "l")
  # matplot(result$state_vecs, type = "l", lty = 2, add = T)

  log_like <- logLik(object = result, data_ = sims$res, id =  sims$res$id)

  # dput(log_like)
  old <- structure(-188.884654279978, class = "logLik")
  expect_equal(log_like, old, tolerance = 1e-3)


  sims <- logit_sim_200
  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ ddFixed(x1) + ddFixed(x2) + . - id - x1 - x2,
    sims$res,
    by = 1,
    control = ddhazard_control(
      n_max = 10^4, eps = 10^-2, est_Q_0 = F, fixed_terms_method = "M_step"),
    Q_0 = diag(1e5, 9),
    Q = diag(1e-2, 9),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F)

  log_like <- logLik(result)

  # matplot(sims$betas, lty = 1, type = "l")
  # matplot(result$state_vecs, type = "l", lty = 2, add = T, col = c(1, 4:11))
  # abline(h = result$fixed_effects, col = 2:3, lty = 2)
  old <- structure(-180.8480151206968,
                   class = "logLik")
  expect_equal(log_like, old, tolerance = 1e-6)

  sims <- exp_sim_200
  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ ddFixed(x1) + ddFixed(x2) + . - id - x1 - x2,
    sims$res,
    by = 1,
    control = ddhazard_control(fixed_terms_method = "E_step"),
    Q_0 = diag(1e5, 9),
    Q = diag(1e-2, 9),
    max_T = 10, model = "exponential",
    id = sims$res$id, order = 1,
    verbose = F)

  log_like <- logLik(result)

  # matplot(sims$betas, lty = 1, type = "l")
  # matplot(result$state_vecs, type = "l", lty = 2, add = T, col = c(1, 4:11))
  # abline(h = result$fixed_effects, col = 2:3, lty = 2)
  # dput(log_like)

  old <- structure(-221.332924990221, class = "logLik")
  expect_equal(log_like, old, tolerance = 1e-6)
})


