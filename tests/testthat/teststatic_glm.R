if(interactive()){
  get_design_matrix <- function(...) environment(ddhazard)$get_design_matrix(...)
}

sims <- test_sim_func_logit(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = .5, re_draw = T, beta_start = 0,
                            intercept_start = -3, sds = c(.1, rep(1, 3)))

test_that("Static glm yields expected number of events, correct rows and same result when supplying risk_obj", {
  form <- survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event
  res <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 10, id = sims$res$id,
    family = "binomial", model = T)

  fail_rows <- sims$res[sims$res$event & sims$res$tstop <= 10, ]
  fail_rows_static <- res$model[res$model$Y == 1, ]

  expect_equal(nrow(fail_rows), nrow(fail_rows_static))
  expect_equal(fail_rows_static$`(weights)`, rep(1, nrow(fail_rows_static)))

  expect_true(any(res$model$`(weights)` > 1))

  design_mat = get_design_matrix(form, sims$res)
  risk_obj = get_risk_obj(Y = design_mat$Y, by = 1, max_T = 10,
                          id = sims$res$id)
  res_own_risk_obj <- dynamichazard::static_glm(
    form = form, data = sims$res, risk_obj = risk_obj)

  expect_equal(res_own_risk_obj$coefficients, res$coefficients)
})

# matplot(sims$betas, type = "l", ylim = range(sims$betas, res$coefficients),
#         col = 1:4)
# sum(sims$res$event)
# abline(h = res$coefficients, col  = 1:4)

set.seed(33587)
sims <- test_sim_func_exp(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
                          x_range = 1, x_mean = .5, re_draw = T, beta_start = 0,
                          intercept_start = -3, sds = c(.1, rep(1, 3)))

test_that("static glm gives results with exponential that match previous computations", {
  form <- survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event
  res <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 10, id = sims$res$id,
    family = "exponential", model = T)

  expect_equal(unname(res$coefficients), c(-3.3734428161003782,  0.6854894359705571,  1.3591615410115157, -0.8936444236277131))

  # test with lower max_T
  res_lower <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 6, id = sims$res$id,
    family = "exponential", model = T)

  expect_equal(unname(res_lower$coefficients), c(-3.04404847708056359, -0.01960283768708104,  0.09531074598781726, -1.53279736136739375))
})

test_that("design_matrix yields equal result with different values of use_weights", {
  form <- formula(survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event, data = sims$res)
  res <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 10, id = sims$res$id,
    family = "logit", model = T)

  expect_equal(unname(res$coefficients), c(-3.0034120780114382, 0.4201450126532649, 0.7693816806898961, -0.5438611125188229))


  data_f <- get_survival_case_Weigths_and_data(formula = form, data = sims$res, by = 1,
                                               max_T = 10, id = sims$res$id, use_weights = F)

  expect_true(all(data_f$weigths == 1))

  form <- update(old = form, new = Y ~ x1 + x2 + x3 , data = data_f)
  glm_res <- glm(form, family = binomial, data = data_f, weights = data_f$weigths)

  expect_equal(glm_res$coefficients, res$coefficients)
})

# cols <- rainbow(4)
# matplot(sims$betas, type = "l", lty = 1, col = cols)
# abline(h = res$coefficients, lty = 2, col = cols)
# sum(sims$res$event)

