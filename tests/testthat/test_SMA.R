context("Testing SMA")

if(interactive()){
  library(survival); library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")

  exp_model_names <- with(environment(ddhazard), exp_model_names)
}

#####
# Test for logit model

test_that("NR method for logit function gives correct values for logit", {
  expect_equal(
    SMA_hepler_logit_compute_length(0, .3, .6, 1, T),
    -0.117687, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_logit_compute_length(0, .3, .6, 1, F),
    -1.34459, tolerance = 1e-5)

  expect_equal(
    SMA_hepler_logit_compute_length(0, .1, .4, 1, T),
    0.222731, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_logit_compute_length(0, .1, .4, 1, F),
    -2.4115, tolerance = 1e-5)

  expect_equal(
    SMA_hepler_logit_compute_length(1, .2, .3, 1, T),
    -0.0518549, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_logit_compute_length(1, .2, .3, 1,  F),
    -1.62284, tolerance = 1e-5)

  expect_equal(
    SMA_hepler_logit_compute_length(1, .3, .1, 5, T),
    0.909057, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_logit_compute_length(1, .3, .1, 5,  F),
    -2.15823, tolerance = 1e-5)
})

cl <- quote(ddhazard(Surv(tstart, tstop, death == 2) ~ age + edema +
               log(albumin) + log(protime) + log(bili), pbc2,
             id = pbc2$id, by = 100, max_T = 3600,
             control = ddhazard_control(
               method = "SMA", posterior_version = "woodbury"),
             Q_0 = diag(rep(100000, 6)), Q = diag(rep(0.01, 6))))

test_that("Logit model for posterior_approx gives previous found values", {
  set.seed(950466)
  f1 <- eval(cl)

  # plot(f1)
  f1 <- f1[c("state_vecs", "state_vecs")]
  expect_known_value(f1, "posterior_approx_logit_pbc.RDS")
})

test_that("Changing between woodbury and cholesky makes a slight difference with PBC",{
  expect_equal(eval(cl$control)$posterior_version, "woodbury")

  set.seed(seed <- 5517547)
  f1 <- eval(cl)
  set.seed(seed)
  f2 <- eval(cl)

  expect_equal(f1[c("state_vars", "state_vecs")], f2[c("state_vars", "state_vecs")],
               tolerance = 1e-16)

  ctrl <- eval(cl$control)
  ctrl$posterior_version <- "cholesky"
  cl$control <- quote(ctrl)
  set.seed(seed)
  f2 <- eval(cl)

  expect_true(is.character(all.equal(
    f1[c("state_vars", "state_vecs")],
    f2[c("state_vars", "state_vecs")], tolerance = 1e-16)))

  # The difference is not that big though
  expect_equal(f1[c("state_vars", "state_vecs")],
               f2[c("state_vars", "state_vecs")], tolerance = 1e-4)
})

test_that("Logit model for posterior_approx differs due to permutation", {
  set.seed(84766)
  f1 <- eval(cl)
  f2 <- eval(cl)

  expect_true(is.character(all.equal(f1$state_vecs, f2$state_vecs, tolerance = 1e-6)))
})

test_that("Logit model for posterior_approx gives previous found values with weights", {
  cl <-  quote(ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(
      est_Q_0 = F, method = "SMA", save_data = F, save_risk_set = F,
      eps = 1e-2),
    Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45, order = 1))

  f1 <- eval(cl) # no weigths

  set.seed(10211)
  ws <- sample(1:3, nrow(head_neck_cancer), replace = T)
  cl$weights <- quote(ws)
  f2 <- eval(cl) # with weigths

  expect_true(!isTRUE(all.equal(
    f1$state_vecs, f2$state_vecs, tolerance = 1e-8)))
  expect_equal(f1$state_vecs, f2$state_vecs, tolerance = 3e-1)
})

test_that("Chaning the learning changes the result for the posterior approx method",{
  skip_on_cran()
  ctrl <- ddhazard_control(est_Q_0 = F, method = "SMA",
                           save_data = F, save_risk_set = F)
  cl <-  quote(ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ctrl,
    Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45, order = 1))

  set.seed(seed <- 685617)
  f1 <- eval(cl)

  set.seed(seed)
  ctrl$LR <- .5
  f2 <- eval(cl)

  f1$control <- NULL
  f1$LR <- NULL
  f1$n_iter <- NULL
  f1$call <- NULL
  f2$control <- NULL
  f2$LR <- NULL
  f2$call <- NULL
  f2$n_iter <- NULL
  expect_true(is.character(all.equal(f1, f2)))
  expect_equal(f1, f2, tolerance = .2)
})

test_that("Second order model gives previous found result for posterior approx", {
  set.seed(1495821)
  f1 <-  ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(method = "SMA"),
    Q_0 = diag(1, 4), Q = diag(0.01, 2),
    max_T = 30, order = 2)

  # plot(f1)
  f1 <- f1[c("state_vecs", "state_vecs")]
  # save_to_test(f1, "posterior_approx_logit_2_order")

  expect_equal(f1, read_to_test("posterior_approx_logit_2_order"), tolerance = 1.490116e-08)
})

test_that("Posterior gives previous found results with large by length for pbc data with logit", {
  set.seed(536705)
  f1 <- ddhazard(Surv(tstart, tstop, death == 2) ~ age + edema +
                  log(albumin) + log(protime) + log(bili), pbc2,
                 id = pbc2$id, by = 300, max_T = 3600,
                 control = ddhazard_control(method = "SMA"),
                 Q_0 = diag(rep(100000, 6)), Q = diag(rep(1e-3, 6)))

  # plot(f1)
  f1 <- f1[c("state_vecs", "state_vecs")]
  # save_to_test(f1, "posterior_approx_logit_pbc_large_by")

  expect_equal(f1, read_to_test("posterior_approx_logit_pbc_large_by"), tolerance = 1e-7)
})

#####
# Test for exponential model

test_that("NR method for logit function gives correct values for Exponential", {
  expect_equal(
    SMA_hepler_exp_compute_length(0, .2, .1, 1, T, 1),
    -0.0733015, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_exp_compute_length(0, .2, .1, 1, F, 1),
    -1.09029, tolerance = 1e-5)

  expect_equal(
    SMA_hepler_exp_compute_length(4, .4, .1, 1, T, 1),
    -2.8445, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_exp_compute_length(4, .4, .1, 1, F, 1),
    -3.12465, tolerance = 1e-5)

  expect_equal(
    SMA_hepler_exp_compute_length(0, 1, .2, 2, T, 1),
    -0.0506302, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_exp_compute_length(0, 1, .2, 2, F, 1),
    -0.631692, tolerance = 1e-5)

  expect_equal(
    SMA_hepler_exp_compute_length(0, .5, .05, 1, T, 10),
    -1.43386, tolerance = 1e-5)
  expect_equal(
    SMA_hepler_exp_compute_length(0, .5, .05, 1, F, 10),
    -1.76385, tolerance = 1e-5)
})


cl <- quote(ddhazard(
  Surv(tstart, tstop, death == 2) ~ age + edema +
    log(albumin) + log(protime) + log(bili), pbc2,
  id = pbc2$id, by = 100, max_T = 3600,
  model = "exp_clip_time_w_jump",
  control = ddhazard_control(method = "SMA", eps = 1e-2),
  Q_0 = diag(rep(100000, 6)), Q = diag(rep(0.001, 6))))

test_that("Exponential model for posterior_approx gives previous found values", {
  set.seed(507958)
  f1 <- eval(cl)

  # plot(f1)
  f1 <- f1[c("state_vecs", "state_vecs")]
  # save_to_test(f1, "posterior_approx_exp_pbc")

  expect_equal(f1, read_to_test("posterior_approx_exp_pbc"), tolerance = 1.490116e-08)
})

test_that("Exponential model yields the same results for all the method inputs with same seed", {
  seed <- 259430

  set.seed(seed)
  cl$model <- "exp_clip_time_w_jump"
  set.seed(seed)
  f1 <- eval(cl)

  for(n in exp_model_names[exp_model_names != "exp_clip_time_w_jump"]){
    set.seed(seed)
    cl$model <- n
    f2 <- eval(cl)

    expect_equal(f1$state_vecs, f2$state_vecs)
  }
})
