context("Testing static_glm")

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
  # print(unname(c(tmp$coefficients)), digits = 16)

  expect_equal(unname(c(tmp$coefficients)),
               c(-3.1239104299217413, 0.5945953199533323, 1.0979243555250271, -0.9534880608516028))

  # test with lower max_T
  res_lower <- dynamichazard::static_glm(
    form = form, data = sims$res, by = 1, max_T = 6, id = sims$res$id,
    family = "exponential", model = T)

  tmp <- res_lower["coefficients"]

  expect_equal(unname(c(tmp$coefficients)),
               c(-2.7762125831981495, -0.3728709821372309, -0.2619269145416432, -1.4856958985013418 ))
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


test_that("get_survival_case_weights_and_data work w/ factors and weights",{
  tmp <- data.frame(
    start = 0:4 * 2,
    stop = 1:5 * 2,
    x = factor(c("a", "a", "b", "b", "c")))

  suppressWarnings(
    res <- get_survival_case_weights_and_data(
      Surv(start, stop, rep(1, 5)) ~ x,
      data = tmp, by = 1, max_T = 10, use_weights = F))

  expect_equal(tmp$x, res$X$x[1:5*2])
})

test_that("Gets error with only_coef = TRUE and no mf",{
  expect_error(
    static_glm(only_coef = T),
    "^mf must be supplied when only_coef = TRUE$")
})

test_that("Gives depreciated warning about speedglm argument", {
  sims <- logit_sim_200$res

  expect_warning(
    static_glm(
      Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
      data = sims, by = 1, max_T = 10, family = "logit",
      id = sims$id, only_coef = FALSE,
      speedglm = TRUE # should give warning
      ),
    ".speedglm. have been depreciated. Use .method_use.")
})

test_that("Gets same with different methods", {
  skip_if_not_installed("speedglm")

  frm <- Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id + x1^2

  for(fam in c("logit", "exponential")){
    sims <- if(fam == "logit")
      logit_sim_200$res else exp_sim_200$res

    get_fit <- function(method_use)
      static_glm(
        frm, data = sims, by = 1, max_T = 10, family = fam,
        id = sims$id, only_coef = FALSE, method_use = method_use)

    f1 <- get_fit("glm")
    f2 <- get_fit("speedglm")

    expect_equal(f1$coefficients, f2$coefficients)

    tmp <- get_design_matrix(frm, sims)
    get_fit <- function(method_use)
      static_glm(
        frm, data = sims, by = 1, max_T = 10, family = fam,
        id = sims$id, only_coef = TRUE, method_use = method_use,
        mf = cbind(tmp$X, tmp$fixed_terms))

    f3 <- get_fit("glm")
    f4 <- get_fit("speedglm")
    f5 <- get_fit("parallelglm")

    expect_equal(f1$coefficients, f3)
    expect_equal(f3, f4, check.attributes = FALSE)
    expect_equal(f4, f5, check.attributes = FALSE)
  }
})
