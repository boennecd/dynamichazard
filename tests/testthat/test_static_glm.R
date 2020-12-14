context("Testing static_glm")

set.seed(615015)
sims <- test_sim_func_logit(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = .5, re_draw = T, beta_start = 0,
                            intercept_start = -3, sds = c(.1, rep(1, 3)))

test_that("Static glm yields expected number of events, correct rows and same result when supplying risk_obj", {
  form <- survival::Surv(tstart, tstop, event) ~ . - id
  res <- static_glm(
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
  res_own_risk_obj <- static_glm(
    form = form, data = sims$res, risk_obj = risk_obj)

  expect_equal(res_own_risk_obj$coefficients, res$coefficients)
})

test_that("static glm gives results with exponential that match previous computations", {
  set.seed(33587)
  sims <- test_sim_func_exp(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = .5, re_draw = T, beta_start = 0,
                            intercept_start = -3, sds = c(.1, rep(1, 3)))

  form <- survival::Surv(tstart, tstop, event) ~ . - id
  res <- static_glm(
    form = form, data = sims$res, by = 1, max_T = 10, id = sims$res$id,
    family = "exponential", model = T)

  tmp <- res["coefficients"]
  # print(unname(c(tmp$coefficients)), digits = 16)

  expect_equal(unname(c(tmp$coefficients)),
               c(-3.1239104299217413, 0.5945953199533323, 1.0979243555250271, -0.9534880608516028))

  # test with lower max_T
  res_lower <- static_glm(
    form = form, data = sims$res, by = 1, max_T = 6, id = sims$res$id,
    family = "exponential", model = T)

  tmp <- res_lower["coefficients"]

  expect_equal(unname(c(tmp$coefficients)),
               c(-2.7762125831981495, -0.3728709821372309, -0.2619269145416432, -1.4856958985013418 ))
})

test_that("design_matrix yields equal result with different values of use_weights", {
  form <- formula(survival::Surv(tstart, tstop, event) ~ . - id, data = sims$res)
  res <- static_glm(
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
    "^mf must be supplied when .only_coef = TRUE.$")
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

  frm <- Surv(tstart, tstop, event) ~ . - id + x1^2

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
    f5 <- get_fit("parallelglm_quick")
    f6 <- get_fit("parallelglm_QR")

    expect_equal(f1$coefficients, f3)
    expect_equal(f3, f4, check.attributes = FALSE)
    expect_equal(f4, f5, check.attributes = FALSE)
    expect_equal(f5, f6, check.attributes = FALSE)
  }
})

test_that("get_survival_case_weights_and_data works when columns names are already used",{
  dat <- data.frame(
    id = c(1, 2, 3),
    Y  = c(TRUE, FALSE, FALSE),
    weights = c(9, 10, 2),
    t  = c(2, 2, 2),
    x  = c(2, 1, 0))

  expect_warning(out <- get_survival_case_weights_and_data(
    Surv(t, Y) ~ x, data = dat, by = 1, max_T = 2, use_weights = FALSE,
    id = dat$id),
    "Column called 'Y' is already in the data.frame. Will use'Y.1'instead", fixed = TRUE)

  expect_equal(
    out$X,
    structure(list(
      id = c(1, 2, 3, 1, 2, 3),
      Y = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE),
      weights = c(9, 10, 2, 9, 10, 2),
      t = c(2, 2, 2, 2, 2, 2),
      x = c(2, 1, 0, 2, 1, 0),
      Y.1 = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), t.1 = c(1, 1, 1, 2, 2, 2),
      weights.1 = c(1, 1, 1, 1, 1, 1)),
      .Names = c("id", "Y", "weights", "t", "x", "Y.1", "t.1", "weights.1"),
      row.names = c(NA, -6L), class = "data.frame"))

  expect_warning(
    fit <- static_glm(Surv(t, Y) ~ -1 + x,
                      data = dat, by = 1, max_T = 2, id = dat$id),
    "Column called 'weights' is already in the data.frame. Will use'weights.1'instead.")

  expect_equal(fit$coefficients,
               glm(Y.1 ~ -1 + x, binomial(), out$X)$coefficients)
})

test_that("get_survival_case_weights_and_data works with observation in counting setup with `holes'",{
  dat <- data.frame(
    id     = c(1 , 1, 2),
    tstart = c(0 , 3, 0),
    tstop  = c(.5, 4, 4),
    event  = c(  F, T, F),
    x      = 1:3)

  res <- get_survival_case_weights_and_data(
    Surv(tstart, tstop, event) ~ x,
    data = dat, by = 1, max_T = 4, id = dat$id,
    use_weights = FALSE)

  expect_equal(
    res$X,
    structure(list(
      id = c(1, 2, 2, 2, 2, 1),
      tstart = c(0, 0, 0, 0, 0, 3),
      tstop = c(0.5, 4, 4, 4, 4, 4),
      event = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
      x = c(1L, 3L, 3L, 3L, 3L, 2L),
      Y = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
      t = c(1, 1, 2, 3, 4, 4),
      weights = c(1, 1, 1, 1, 1, 1)),
      .Names = c("id", "tstart", "tstop", "event", "x", "Y", "t", "weights"),
      row.names = c(NA, -6L), class = "data.frame"))

  res <- get_survival_case_weights_and_data(
    Surv(tstart, tstop, event) ~ x,
    data = dat, by = 1, max_T = 4, id = dat$id,
    use_weights = TRUE)

  expect_equal(
    res$X,
    structure(list(
      Y = c(0, 0, 1), id = c(1, 2, 1),
      tstart = c(0, 0, 3), tstop = c(0.5, 4, 4),
      event = c(FALSE, FALSE, TRUE), x = c(1L, 3L, 2L), weights = c(1, 4, 1)),
      .Names = c("Y", "id", "tstart", "tstop", "event", "x", "weights"),
      row.names = c(NA, -3L), class = "data.frame"))

  dat <- data.frame(
    id     = c(1, 1  , 1  , 2),
    tstart = c(0, 1.5, 3.5, 0),
    tstop  = c(1, 2  , 4  , 4),
    event  = c(F, F  , T  , F),
    x      = 1:4)

  res <- get_survival_case_weights_and_data(
    Surv(tstart, tstop, event) ~ x,
    data = dat, by = 1, max_T = 4, id = dat$id,
    use_weights = FALSE)

  expect_equal(
    res$X,
    structure(list(
      id = c(1, 2, 2, 2, 2), tstart = c(0, 0, 0, 0, 0),
      tstop = c(1, 4, 4, 4, 4), event = c(FALSE, FALSE, FALSE, FALSE, FALSE),
      x = c(1L, 4L, 4L, 4L, 4L), Y = c(FALSE, FALSE, FALSE, FALSE, FALSE),
      t = c(1, 1, 2, 3, 4), weights = c(1, 1, 1, 1, 1)),
      .Names = c("id", "tstart", "tstop", "event", "x", "Y", "t", "weights"),
      row.names = c(NA, -5L), class = "data.frame"))
})
