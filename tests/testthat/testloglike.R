if(interactive())
  library(testthat); library(survival); library(parallel); source("R/test_utils.R")

test_that("Verbose on ddhazard prints a log likelihood", {
  expect_output({
    result = ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1,
      control = list(n_max = 10^4, eps = 10^-4, est_Q_0 = F),
      a_0 = rep(0, 2), Q_0 = diag(1, 2),
      max_T = 45,
      id = head_neck_cancer$id, order = 1,
      verbose = 5
    )
  }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood is\\s+")

  expect_output({
    result = ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1,
      control = list(n_max = 10^4, eps = 10^-4, est_Q_0 = F),
      a_0 = rep(0, 2), Q_0 = diag(1, 2),
      max_T = 45,
      id = head_neck_cancer$id, order = 1,
      verbose = 5,
      model = "exponential"
    )
  }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood is\\s+")

  expect_output({
    result = ddhazard(
      formula = survival::Surv(start, stop, event) ~ ddFixed(group),
      data = head_neck_cancer,
      by = 1,
      control = list(n_max = 10^4, eps = 10^-4, est_Q_0 = F),
      a_0 = 0, Q_0 = as.matrix(1),
      max_T = 45,
      id = head_neck_cancer$id, order = 1,
      verbose = 5,
      model = "exponential"
    )
  }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood is\\s+")
})

test_that("logLik for head_neck_cancer data set match previous results", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(n_max = 10^4, eps = 10^-4, est_Q_0 = F),
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1,
    verbose = F
  )

  log_like <- logLik(object = result)
  result$data <- NULL
  result$risk_set <- NULL
  logLik(object = result, data = head_neck_cancer, id = head_neck_cancer$id)

  old <- structure(-340.3724737164279,
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
    a_0 = rep(0, 4), Q_0 = diag(c(1, 1, 1, 1)),
    Q = diag(c(1e-4, 1e-4, 0, 0)),
    control = list(est_Q_0 = F, n_max = 10^4, eps = 10^-4),
    max_T = 45,
    id = head_neck_cancer$id, order = 2,
    verbose = F
  )

  log_like <- logLik(object = result, data_ = head_neck_cancer, id =  head_neck_cancer$id)

  old <- structure(-349.569928801165,
                   class = "logLik",
                   df = 2 * 2 + 3)

  expect_equal(log_like, old)
})

test_that("logLik for head_neck_cancer data set match previous results with fixed effects", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1,
    control = list(n_max = 10^4, eps = 10^-4, est_Q_0 = F),
    a_0 = 0, Q_0 = as.matrix(1),
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  log_like <- logLik(result, data = head_neck_cancer)

  old <- structure(-302.8375676988051,
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
    id = head_neck_cancer$id, order = 1))

  tmp_design <- get_survival_case_Weigths_and_data(
    formula = form, data = head_neck_cancer, by = 1, max_T = 45, id = head_neck_cancer$id,
    use_weights = F)

  glm_fit <- glm(Y ~ as.factor(group), binomial(), tmp_design)

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
    by = 1,
    control = list(n_max = 10^4, eps = 10^-2, est_Q_0 = F),
    a_0 = rep(0, 5), Q_0 = diag(1, 5),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F
  )

  log_like <- logLik(object = result)
  old <- structure(-2414.852657205478,
                   class = "logLik",
                   df = 5 + 5 * (1 + 5) / 2)
  expect_equal(log_like, old)


  set.seed(35374)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 5, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -3, sds = c(.1, rep(1, 5)))

  result = ddhazard(
    survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    sims$res,
    by = 1,
    control = list(n_max = 10^4, eps = 10^-2, est_Q_0 = F),
    a_0 = rep(0, 6), Q_0 = diag(10, 6),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F, model = "exponential"
  )

  log_like <- logLik(object = result, data_ = sims$res, id =  sims$res$id)
  old <- structure(-1561.570688141223,
                   class = "logLik",
                   df = 6 + 6 * (1 + 6) / 2)
  expect_equal(log_like, old, tolerance = 1e-6)





  set.seed(4578234)
  sims <- test_sim_func_logit(n_series = 2e3, n_vars = 5, t_0 = 0, t_max = 10,
                              x_range = 1, x_mean = -.2, re_draw = T)

  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ ddFixed(x1) + ddFixed(x2) + x3 + x4 + x5,
    sims$res,
    by = 1,
    control = list(n_max = 10^4, eps = 10^-2, est_Q_0 = F),
    a_0 = rep(0, 4), Q_0 = diag(1, 4),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F)

  log_like <- logLik(result)
  old <- structure(-2909.663034719353,
                   class = "logLik",
                   df = 4 + 4 * (1 + 4) / 2 + 2)
  expect_equal(log_like, old, tolerance = 1e-6)
})

test_that("Add one example of exponential fit and loglike with fixed effects", expect_true(F))


