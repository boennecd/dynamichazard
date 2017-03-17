if(interactive()){
  library(dynamichazard)
  library(testthat)
  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")
  get_design_matrix <- function(...) environment(ddhazard)$get_design_matrix(...)
}


# Had issues with win builder. Thus, these lines
test_name <- "static_glm"
cat("\nRunning", test_name, "\n")



set.seed(615015)
sims <- test_sim_func_logit(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
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

test_that("static glm gives results with exponential that match previous computations", {
  set.seed(33587)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = .5, re_draw = T, beta_start = 0,
                            intercept_start = -3, sds = c(.1, rep(1, 3)))

  form <- survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event
  res <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 10, id = sims$res$id,
    family = "exponential", model = T)

  tmp <- res["coefficients"]

  expect_equal(unname(c(tmp$coefficients)),
               c(-3.0654969893646609, 0.4515918710312326, 1.0892116143756236, -0.7820465199455019 ))

  # test with lower max_T
  res_lower <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 6, id = sims$res$id,
    family = "exponential", model = T)

  tmp <- res_lower["coefficients"]
  # get_expect_equal(tmp)

  expect_equal(unname(c(tmp$coefficients)),
               c(-2.7659195344671597, -0.6083779316590273, -0.1280512349287461, -1.0481131739738363 ))
})

test_that("design_matrix yields equal result with different values of use_weights", {
  form <- formula(survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event, data = sims$res)
  res <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 10, id = sims$res$id,
    family = "logit", model = T)

  res <- res[c("coefficients")]
  # save_to_test(res, "static1")

  expect_equal(res, read_to_test("static1"))

  data_f <- get_survival_case_weights_and_data(formula = form, data = sims$res, by = 1,
                                               max_T = 10, id = sims$res$id, use_weights = F)$X

  expect_true(all(data_f$weigths == 1))

  form <- update(old = form, new = Y ~ x1 + x2 + x3 , data = data_f)
  glm_res <- glm(form, family = binomial, data = data_f, weights = data_f$weigths)

  expect_equal(glm_res$coefficients, res$coefficients)
})



# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
