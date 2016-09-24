test_that("Verbose on ddhazard prints a log likelihood", {
  expect_output({
    result = ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1,
      n_max = 10^4, eps = 10^-4,
      a_0 = rep(0, 2), Q_0 = diag(1, 2),
      est_Q_0 = F,
      max_T = 45,
      id = head_neck_cancer$id, order_ = 1,
      verbose = 5
    )
  }, regexp = "Iteration\\s+\\d+\\sended with conv criteria\\s+\\d+.\\d+\\s+The log likelihood is\\s+")

  expect_output({
    result = ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1,
      n_max = 10^4, eps = 10^-4,
      a_0 = rep(0, 2), Q_0 = diag(1, 2),
      est_Q_0 = F,
      max_T = 45,
      id = head_neck_cancer$id, order_ = 1,
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
    n_max = 10^4, eps = 10^-4,
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    est_Q_0 = F,
    max_T = 45,
    id = head_neck_cancer$id, order_ = 1,
    verbose = F
  )

  log_like <- logLik(object = result, data_ = head_neck_cancer, id =  head_neck_cancer$id)

  old <- structure(-340.3724737164279,
                   class = "logLik",
                   df = 2 + 3)

  expect_equal(log_like, old)
})

test_that("logLik for head_neck_cancer data set with second order model", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    n_max = 10^4, eps = 10^-4,
    a_0 = rep(0, 4), Q_0 = diag(c(1, 1, 1, 1)),
    Q = diag(c(1e-4, 1e-4, 0, 0)),
    est_Q_0 = F,
    max_T = 45,
    id = head_neck_cancer$id, order_ = 2,
    verbose = F
  )

  log_like <- logLik(object = result, data_ = head_neck_cancer, id =  head_neck_cancer$id)

  old <- structure(-349.569928801165,
                   class = "logLik",
                   df = 2 * 2 + 3)

  expect_equal(log_like, old)
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
    n_max = 10^4, eps = 10^-2,
    a_0 = rep(0, 5), Q_0 = diag(1, 5),
    est_Q_0 = F,
    max_T = 10,
    id = sims$res$id, order_ = 1,
    verbose = F
  )

  log_like <- logLik(object = result, data_ = sims$res, id =  sims$res$id)
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
    n_max = 10^4, eps = 10^-2,
    a_0 = rep(0, 6), Q_0 = diag(10, 6),
    est_Q_0 = F,
    max_T = 10,
    id = sims$res$id, order_ = 1,
    verbose = F, model = "exponential"
  )

  log_like <- logLik(object = result, data_ = sims$res, id =  sims$res$id)
  old <- structure(-1561.570688141223,
                   class = "logLik",
                   df = 6 + 6 * (1 + 6) / 2)
  expect_equal(log_like, old, tolerance = 1e-6)
})


