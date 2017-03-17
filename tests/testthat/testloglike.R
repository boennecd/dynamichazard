if(interactive()){
  library(testthat); library(survival); library(parallel)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")

  library(dynamichazard)
  exp_model_names <- with(environment(ddhazard), exp_model_names)
}


# Had issues with win builder. Thus, these lines
test_name <- "loglike"
cat("\nRunning", test_name, "\n")

test_that("verbose prints log likelihood",{
  for(method in c("EKF", "UKF")){
    for(o in c(1, 2)){
      Q_0_arg <- if(o == 1) diag(1, 2) else diag(c(1, 1, .1, .1))
      for(m in c(exp_model_names, "logit")){
        if(m == "exp_clip_time" && method == "UKF")
          next
        expect_output({
          ddhazard(survival::Surv(start, stop, event) ~ group,
                   data = head_neck_cancer, id = head_neck_cancer$id,
                   by = 2, max_T = 30, Q = diag(.01, 2),
                   Q_0 = Q_0_arg, model = m, order = o,
                   verbose = 5,
                   control = list (eps = .1, method = method))
        }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood is\\s+")
}}}})

test_that("verbose prints log likelihood with fixed effects",{
  for(method in c("EKF", "UKF")){
    for(fixed_method in c("M_step", "E_step")){
      Q_0_arg <- if(method == "EKF") 1e5 else 1
      for(m in c(exp_model_names, "logit")){
        if(m == "exp_clip_time" && method == "UKF")
          next
        expect_output({
          ddhazard(survival::Surv(start, stop, event) ~ ddFixed(group),
                   data = head_neck_cancer, id = head_neck_cancer$id,
                   by = 2, max_T = 30, Q = diag(.1, 1),
                   Q_0 = Q_0_arg, model = m,
                   verbose = 5,
                   control = list (eps = .1, method = method,
                                   fixed_terms_method = fixed_method))
        }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood is\\s+")
      }}}})

test_that("logLik for head_neck_cancer data set match previous results", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F),
    a_0 = rep(0, 2), Q_0 = diag(1e5, 2),
    Q = diag(1e-2, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1,
    verbose = F)

  log_like <- logLik(object = result)
  result$data <- NULL
  result$risk_set <- NULL
  logLik(object = result, data = head_neck_cancer, id = head_neck_cancer$id)

  # plot(result)
  # print(log_like, digits = 16)

  old <- structure(-149.0302744634669,
                   class = "logLik",
                   df = 2 + 3)

  expect_equal(log_like, old)
})

test_that("Saving or not saving risk set or data gives the same result", {
  arg_list <- list(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1,
    verbose = F)

  control_fit <- do.call(ddhazard, arg_list)
  control_fit$logLik <- logLik(control_fit)


  for(save_risk_set in c(T, F))
    for(save_data in c(T, F)){
      if(save_risk_set && save_data)
        next

      arg_list$control <- list(save_risk_set = save_risk_set, save_data = save_data)
      new_fit <- do.call(ddhazard, arg_list)
      new_fit$logLik <- logLik(new_fit, data = if(save_data) NULL else head_neck_cancer,
                               id = if(save_risk_set) NULL else head_neck_cancer$id)

      expect_equal(new_fit$logLik, control_fit$logLik)
    }
})

test_that("logLik for head_neck_cancer data set with second order model", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    a_0 = rep(0, 4), Q_0 = diag(1000000, 4),
    Q = diag(1e-2, 2),
    control = list(est_Q_0 = F, n_max = 10^4, eps = 10^-2),
    max_T = 45,
    id = head_neck_cancer$id, order = 2,
    verbose = F
  )

  log_like <- logLik(object = result, data_ = head_neck_cancer, id =  head_neck_cancer$id)


  # plot(result)
  # print(log_like, digits = 16)

  old <- structure(-129.1952699930929,
                   class = "logLik",
                   df = 2 * 2 + 3)

  expect_equal(log_like, old)
})

test_that("logLik for head_neck_cancer data set match previous results with fixed effects", {
  suppressWarnings(result <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, LR = .8,
                   fixed_terms_method = "M_step"),
    a_0 = 0, Q_0 = 1e5, Q = 1e-1,
    max_T = 45,
    id = head_neck_cancer$id, order = 1))

  # plot(result)
  log_like <- logLik(result, data = head_neck_cancer)

  # print(log_like, digits = 16)

  old <- structure(-247.0969775644861,
                   class = "logLik",
                   df = 1 + 1 + 1)

  expect_equal(log_like, old)
})

test_that("logLik for head_neck_cancer data with only fixed match bigglm", {
  form <- survival::Surv(start, stop, event) ~ -1 + ddFixed(rep(1, length(group))) +
    ddFixed(as.numeric(group == 1))

  suppressWarnings(result <- ddhazard(
    formula = form,
    data = head_neck_cancer,
    by = 1,
    max_T = 45,
    id = head_neck_cancer$id, order = 1,
    control = list(fixed_terms_method = "M_step")))

  tmp_design <- get_survival_case_weights_and_data(
    formula = form, data = head_neck_cancer, by = 1, max_T = 45, id = head_neck_cancer$id,
    use_weights = F)

  glm_fit <- glm(Y ~ as.factor(group), binomial(), tmp_design$X)

  tmp <- logLik(glm_fit)
  attributes(tmp)
  attr(tmp, "nobs") <- NULL

  expect_equal(c(unname(result$fixed_effects)), unname(glm_fit$coefficients))
  expect_equal(logLik(result), tmp)
})


##############

test_that("logLik for simulated data versus old results", {
  set.seed(21280)
  sims <- test_sim_func_logit(n_series = 2e3, n_vars = 5, t_0 = 0, t_max = 10,
                              x_range = 1, x_mean = -.2, re_draw = T)
  sims$res <- as.data.frame(sims$res)

  result = ddhazard(
    survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event - 1,
    sims$res,
    by = 1, Q = diag(1e-1, 5),
    a_0 = rep(0, 5), Q_0 = diag(1e5, 5),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F
  )

  # matplot(sims$betas, lty = 1, type = "l")
  # matplot(result$state_vecs, type = "l", lty = 2, add = T)

  log_like <- logLik(object = result)

  # print(log_like, digits = 16)

  old <- structure(-2563.452319960646,
                   class = "logLik",
                   df = 5 + 5 * (1 + 5) / 2)
  expect_equal(log_like, old)


  set.seed(35374)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 5, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -3, sds = c(.1, rep(1, 5)))

  suppressMessages(result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    sims$res,
    by = 1, Q = diag(1e-1, 6),
    a_0 = rep(0, 6), Q_0 = diag(1e5, 6),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F, model = "exp_clip_time_w_jump"))

  # matplot(sims$betas, lty = 1, type = "l")
  # matplot(result$state_vecs, type = "l", lty = 2, add = T)

  log_like <- logLik(object = result, data_ = sims$res, id =  sims$res$id)

  # print(log_like, digits = 16)

  old <- structure(-1511.394013234687,
                   class = "logLik",
                   df = 6 + 6 * (1 + 6) / 2)
  expect_equal(log_like, old, tolerance = 1e-3)





  set.seed(4578234)
  sims <- test_sim_func_logit(n_series = 2e3, n_vars = 5, t_0 = 0, t_max = 10,
                              x_range = 1, x_mean = -.2, re_draw = T)

  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ ddFixed(x1) + ddFixed(x2) + x3 + x4 + x5,
    sims$res,
    by = 1,
    control = list(n_max = 10^4, eps = 10^-2, est_Q_0 = F,
                   fixed_terms_method = "M_step"),
    a_0 = rep(0, 4), Q_0 = diag(1e5, 4),
    Q = diag(1e-1, 4),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F)



  log_like <- logLik(result)

  # print(log_like, digits = 16)

  old <- structure(-2897.877470351145,
                   class = "logLik",
                   df = 4 + 4 * (1 + 4) / 2 + 2)
  expect_equal(log_like, old, tolerance = 1e-6)




  set.seed(4578234)
  sims <- test_sim_func_exp(n_series = 5e2, n_vars = 5, t_0 = 0, t_max = 10,
                            beta_start = (seq_len(5) - 4) / 2, intercept_start = -3,
                            x_range = 2, x_mean = 0, re_draw = T, is_fixed = c(1:2))

  suppressMessages(result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ ddFixed(1) + ddFixed(x1) + x2 + x3 + x4 + x5,
    sims$res,
    by = 1,
    control = list(fixed_terms_method = "E_step"),
    Q_0 = diag(1e5, 4),
    Q = diag(1e-2, 4),
    max_T = 10, model = "exp_clip_time_w_jump",
    id = sims$res$id, order = 1,
    verbose = F))

  log_like <- logLik(result)

  # print(log_like, digits = 16)
  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(result$state_vecs, type = "l", lty = 2, col = 3:6, add = T)
  # abline(h = result$fixed_effects, lty = 2, col = 1:2)

  old <- structure(-509.6158134865298 ,
                   class = "logLik",
                   df = 4 + 4 * (1 + 4) / 2 + 2)
  expect_equal(log_like, old, tolerance = 1e-6)
})


# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")


