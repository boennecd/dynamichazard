context("Testing ddhazard w/ generic things and w/ the the EKF")

# Test on data set that is one of Farhmiers papers
result <- ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  control = ddhazard_control(
    est_Q_0 = F, save_data = F, save_risk_set = F),
  a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
  max_T = 45,
  id = head_neck_cancer$id, order = 1)

result <- result[c("state_vecs", "state_vars","lag_one_cov", "model")]

test_that("get previous results with head_neck", {
  expect_equal(result, read_to_test("ddhazard_head_neck"))
})

test_that("Invalid penalty terms throw error", {
  expect_error(
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, # Use by month intervals
      control = ddhazard_control(denom_term = 0)),
    regexp = "Method not implemented with penalty term 'control\\$denom_term' equal to 0")

  expect_error(
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, # Use by month intervals
      control = ddhazard_control(denom_term = -1)),
    regexp = "Method not implemented with penalty term 'control\\$denom_term' equal to -1")
})

test_that("Get expected warning when no Q or Q_0 is passed", {
  cl <- quote(ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(
      est_Q_0 = F, save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1))

  tmp <- cl
  tmp$Q_0 <- NULL
  expect_warning(
    eval(tmp),
    paste(sQuote("Q_0"), "not supplied. It has been set to a diagonal matrix with diagonal entries equal to"))

  tmp <- cl
  tmp$Q <- NULL
  expect_warning(
    eval(tmp),
    paste(sQuote("Q"), "not supplied. It has been set to a diagonal matrix with diagonal entries equal to"))
})



test_that("Changing convergence criteria change output", {
  ctrl <- ddhazard_control(criteria = "delta_coef", eps = .002)
  cl <- quote(ddhazard(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, # Use by month intervals
    id = head_neck_cancer$id,
    Q_0 = diag(1e5, 2), Q = diag(.1, 2),
    control = ctrl))

  suppressMessages(res1 <- eval(cl))
  ctrl$criteria <- "delta_likeli"
  suppressMessages(res2 <- eval(cl))

  expect_true(res1$n_iter != res2$n_iter)
})

result_exp <- ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1, Q_0 = diag(10000, 2),
  Q = diag(1e-3, 2), a_0 = c(0, 0),
  max_T = 30,
  id = head_neck_cancer$id, order = 1,
  model = "exponential")

test_that("exponential model and logit moels hazzard functions differs", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(
      est_Q_0 = F, save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  expect_true(result_exp$model != result$model)
  expect_true(result_exp$family$name() != result$family$name())
})

test_that("Testing names of output from ddhazard on head and neck cancer dataset", {
  expect_equal(colnames(result_exp$state_vecs), c("(Intercept)", "group1"))
  expect_equal(unlist(dimnames(result_exp$state_vars)), unlist(list(c("(Intercept)", "group1"), c("(Intercept)", "group1"), NULL)))
  expect_equal(unlist(dimnames(result_exp$Q)), rep(c("(Intercept)", "group1"), 2))
  expect_equal(unlist(dimnames(result_exp$Q_0)), rep(c("(Intercept)", "group1"), 2))
})

# plot(result_exp)
result_exp <- result_exp[
  names(result_exp) %in%
    c("state_vecs", "state_vars", "fixed_effects", "Q")]
test_that("Result of exponential model on head_neck_data match previous results", {
  expect_known_value(
    result_exp, file = "head_neck_exp.RDS", update = FALSE,
    check.attributes = FALSE)
})
rm(result_exp)


test_that("Unmacthed control variable throw error",
          expect_error({
            result = ddhazard(
              formula = survival::Surv(start, stop, event) ~ group,
              data = head_neck_cancer,
              by = 1, # Use by month intervals
              a_0 = rep(0, 2), Q_0 = diag(1, 2), # Initial value
              max_T = 45,
              id = head_neck_cancer$id, order = 1,
              control = ddhazard_control(None_existing_parem = 1)
            )}, regexp = "Unused arguments passed to 'ddhazard_control'"))

test_that("Various ways of passing control gives the same but some with warnings", {
  cl <- quote(ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, save_data = F),
    a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1))

  expect_warning(
    f1 <- eval(cl),
    "'ddhazard_control' instead of 'list' to the 'control' argument")

  ctrl <- eval(cl$control)
  cl$control <- quote(ctrl)
  expect_silent(f2 <- eval(cl))
  expect_equal(f1[names(f1) != "call"], f2[names(f2) != "call"])

  cl$control <- quote(ddhazard_control(est_Q_0 = F, save_data = F))
  expect_silent(f3 <- eval(cl))
  expect_equal(f1[names(f1) != "call"], f3[names(f3) != "call"])

  cl$control <- ddhazard_control(est_Q_0 = F, save_data = F)
  expect_silent(f4 <- eval(cl))
  expect_equal(f1[names(f1) != "call"], f4[names(f4) != "call"])
})

test_that("Different non-integer time_scales gives the same result with ddhazard", {
  skip_on_cran()
  scales <- exp(seq(-2, 2, .05))

  fit_exp <- expression(
    ddhazard(
      formula = survival::Surv(start * .by, stop * .by, event) ~ group,
      data = head_neck_cancer,
      by = .by,
      control = ddhazard_control(est_Q_0 = F, save_data = F),
      a_0 = rep(0, 2), Q_0 = diag(1e2, 2), Q = diag(1e-2 / .by, 2),
      id = head_neck_cancer$id, order = 1))

  .by <- scales[1]
  f1 <- eval(fit_exp)
  for(.by in scales[-1]){
    f2 <- eval(fit_exp)
    info <- paste("by =", .by)
    expect_equal(f1$risk_set$risk_sets, f2$risk_set$risk_sets, info = info)
    expect_equal(f1$risk_set$is_event_in, f2$risk_set$is_event_in, info = info)
    expect_equal(f1$state_vecs, f2$state_vecs, info = info)
  }
})

test_that(
  "Old expoential models gives the same results and yields expected message", {
    cl <- quote(ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(10000, 2),
      Q = diag(1e-3, 2), a_0 = c(0, 0),
      max_T = 30,
      id = head_neck_cancer$id, order = 1,
      control = ddhazard_control(eps = .1)))

    cl$model <- "exponential"
    expect_silent(f1 <- eval(cl))

    for(m in c("exp_bin", "exp_clip_time", "exp_clip_time_w_jump")){
      eval(bquote({
        cl$model <- .(m)
        expect_message(
          f2 <- eval(cl),
          .(paste0(".", m, ". is not used after version 0\\.5\\.0\\.")))
        expect_equal(f1[c("state_vars", "state_vecs")],
                     f2[c("state_vars", "state_vecs")],
                     info = .(m))
      }))
    }

  })

test_that("est_a_0 fixes the time zero value", {
  # TODO: make a smarter way to test this...
  sink(tmp_file <- tempfile())
  tryCatch({
    res <- ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(10000, 2),
      Q = diag(1e-3, 2), a_0 = c(0, 0),
      max_T = 30, control = ddhazard_control(est_a_0 = FALSE, debug = TRUE),
      id = head_neck_cancer$id, order = 1,
      model = "logit")

    log_f <- paste0(readLines(tmp_file), collapse = "\n")
    expect_true(grepl(
      "it\\s+1,\\ Starting EM:\\ a_0\\n[^\\n]+\\s+0\\s+0", log_f, perl = TRUE))
    expect_true(grepl(
      paste0("it\\s+", res$n_iter - 1L,
             ",\\ Starting EM:\\ a_0\\n[^\\n]+\\s+0\\s+0"), log_f,
      perl = TRUE))
  }, finally = {
    sink()
    unlink(tmp_file)
  })
})

########
# Test on simulated data

test_that("Result of exponential model gives previous results w/ simulated data", {
  result_exp <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id,
    data = exp_sim_500$res,
    by = 1,
    Q_0 = diag(1e5, 11),
    Q = diag(1e-3, 11),
    control = ddhazard_control(
      save_data = F, save_risk_set = F,
      method = "EKF"),
    max_T = 10,
    id = exp_sim_500$res$id, order = 1,
    model = "exponential")

  # matplot(exp_sim_500$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, lty = 2, type = "l", add = T)
  result_exp <- result_exp[c("state_vars", "state_vecs", "Q")]
  expect_known_value(
    result_exp, "sim_exp.RDS", tolerance = 1e-5, update = FALSE)
})


test_that("Permutating data does not change the results", {
  cl <- quote(ddhazard(
    Surv(stop, event) ~ group, head_neck_cancer,
    by = 1, max_T = 40,
    Q_0 = diag(rep(10000, 2)), Q = diag(rep(0.1, 2))))

  r1 <- eval(cl)
  cl$control <- quote(ddhazard_control(permu = T))
  r2 <- eval(cl)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(!isTRUE(all.equal(
    r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))

  #####
  # With fixed effects
  cl <- quote(ddhazard(
    Surv(tstart, tstop, death == 2) ~ age + ddFixed(edema) +
      log(albumin) + log(protime) + log(bili), pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    Q_0 = diag(rep(10000, 5)), Q = diag(rep(0.001, 5)),
    control = ddhazard_control(n_threads = 1)))

  r1 <- eval(cl)
  cl$control <- quote(ddhazard_control(n_threads = 1, permu = TRUE))
  r2 <- eval(cl)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(!isTRUE(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))

  #####
  # With weigths
  set.seed(94884214)
  w <- sample(1:3, nrow(pbc2), replace = T)

  cl <- quote(ddhazard(
    Surv(tstart, tstop, death == 2) ~ age + edema +
      log(albumin) + log(protime) + log(bili), pbc2,
    id = pbc2$id, by = 100, max_T = 3000,
    weights = w,
    Q_0 = diag(rep(10000, 6)), Q = diag(rep(0.001, 6)),
    control = ddhazard_control(n_threads = 1, LR = .6)))

  r1 <- eval(cl)
  cl$control <- quote(ddhazard_control(n_threads = 1, LR = .6, permu = TRUE))
  r2 <- eval(cl)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(is.character(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))
})

