context("Running tests for GMA")

test_that("Gives same results w/ 1. order logit model", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(method = "GMA"),
    Q_0 = diag(1, 2), Q = diag(0.01, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  # plot(result)
  result <- result[c("state_vecs", "state_vars")]
  # save_to_test(result, "GMA_head_neck_cancer_o1")
  expect_equal(result, read_to_test("GMA_head_neck_cancer_o1"), tolerance = 1.490116e-08)
})

test_that("GMA gives the same w/ all exponential model inputs", {
  args <- list(
    formula = Surv(tstart, tstop, death == 2) ~ edema + log(albumin) + log(bili) + log(protime),
    data = pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = list(method = "GMA"),
    Q_0 = diag(rep(1, 5)), Q = diag(rep(1e-4, 5)))

  args$model <- "exp_bin"

  f1 <- do.call(ddhazard, args)

  # plot(f1)
  f1 <- f1[c("state_vecs", "state_vars")]
  # save_to_test(f1, "GMA_pbc_o1_exp")
  expect_equal(f1, read_to_test("GMA_pbc_o1_exp"), tolerance = 1.490116e-08)

  for(n in exp_model_names[exp_model_names != "exp_bin"]){
    args$model <- n
    f2 <- do.call(ddhazard, args)
    expect_equal(f1$state_vecs, f2$state_vecs)
    expect_equal(f1$state_vars, f2$state_vars)
  }
})

test_that("GMA works w/ second order random walk", {
  f1 <- ddhazard(
    formula = Surv(tstart, tstop, death == 2) ~ edema + log(albumin) + log(bili) + log(protime),
    data = pbc2,
    order = 2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = list(method = "GMA", LR = .33, GMA_NR_eps = 1e-3),
    Q_0 = diag(rep(1, 10)), Q = diag(rep(1e-4, 5)))

  # plot(f1)
  f1 <- f1[c("state_vecs", "state_vars", "Q")]
  # save_to_test(f1, "GMA_pbc_o2_logit")
  expect_equal(f1, read_to_test("GMA_pbc_o2_logit"), tolerance = 1.490116e-08)
})

test_that("GMA works w/ fixed effects in E-step", {
  fit <- ddhazard(
    formula = Surv(tstart, tstop, death == 2) ~ edema + ddFixed(age) + log(albumin) + log(bili) + log(protime),
    data = pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = list(method = "GMA",
                   fixed_terms_method = 'E_step'),
    Q_0 = diag(rep(1, 5)), Q = diag(rep(1e-4, 5)))

  # plot(fit)
  # fit$fixed_effects
  fit <- fit[c("state_vecs", "state_vars", "fixed_effects", "Q")]
  # save_to_test(fit, "GMA_pbc_fixed_E_step")
  expect_equal(fit, read_to_test("GMA_pbc_fixed_E_step"), tolerance = 1.490116e-08)
})

test_that("GMA works w/ fixed effects in M-step", {
  fit <- ddhazard(
    formula = Surv(tstart, tstop, death == 2) ~
      edema + ddFixed(age) + log(albumin) + log(bili) + log(protime),
    data = pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = list(method = "GMA", eps = 1e-3, criteria = "delta_likeli",
                   fixed_terms_method = 'M_step'),
    Q_0 = diag(rep(1, 5)), Q = diag(1e-4, 5))

  # plot(fit)
  # fit$fixed_effects
  fit <- fit[c("state_vecs", "state_vars", "fixed_effects")]
  # save_to_test(fit, "GMA_pbc_fixed_M_step")
  expect_equal(fit, read_to_test("GMA_pbc_fixed_M_step"), tolerance = 1.490116e-08)
})

test_that("Changing hyper parameters w/ GMA changes the result", {
  args <- list(
    formula = Surv(tstart, tstop, death == 2) ~ edema + age + log(albumin) + log(bili) + log(protime),
    data = pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = list(method = "GMA"),
    Q_0 = diag(rep(1, 6)), Q = diag(rep(1e-4, 6)))

  args$control$GMA_NR_eps <- 1
  f1 <- do.call(ddhazard, args)
  args$control$GMA_NR_eps <- .01
  f2 <- do.call(ddhazard, args)

  expect_true(is.character(all.equal(f1$state_vecs, f2$state_vars)))
})

test_that("GMA makes one mesasge when global scoring did not succed within given number of iterations", {
  args <- list(
    formula = Surv(tstart, tstop, death == 2) ~ edema + age + log(albumin) + log(bili) + log(protime),
    data = pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = list(method = "GMA", LR = .5),
    Q_0 = diag(rep(1, 6)), Q = diag(rep(1e-4, 6)))

  args$control$GMA_max_rep <- 1
  expect_warning(
    f1 <- do.call(ddhazard, args),
    "^Failed once to make correction step")

  args$control$GMA_max_rep <- 25
  expect_silent(suppressMessages(f1 <- do.call(ddhazard, args)))
})

test_that("Changing learning rates changes the result w/ GMA", {
  args <- list(
    formula = Surv(tstart, tstop, death == 2) ~ edema + age + log(albumin) + log(bili) + log(protime),
    data = pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    control = list(method = "GMA", LR =.75),
    Q_0 = diag(rep(1, 6)), Q = diag(rep(1e-4, 6)))

  f1 <- do.call(ddhazard, args)
  args$control$LR <- .5
  f2 <- do.call(ddhazard, args)

  expect_true(is.character(all.equal(f1$state_vecs, f2$state_vecs)))
})

test_that("GAM works w/ weights", {
  args <-  list(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, method = "GMA"),
    Q_0 = diag(1, 2), Q = diag(0.01, 2),
    max_T = 45, order = 1)

  f1 <- do.call(ddhazard, args) # no weigths

  set.seed(10211)
  ws <- sample(1:3, nrow(head_neck_cancer), replace = T)
  args$weights = ws
  f2 <- do.call(ddhazard, args) # with weigths

  expect_true(is.character(all.equal(
    f1$state_vecs, f2$state_vecs, tolerance = 1e-8)))
  expect_equal(
    f1$state_vecs, f2$state_vecs, tolerance = 1)

  dum_dat <-
    head_neck_cancer[unlist(mapply(rep, x = seq_along(ws), times = ws)), ]
  args$weights <- NULL
  args$data <- dum_dat
  f3 <- do.call(ddhazard, args) # with dummy data mimic weigths

  expect_equal(f2$state_vecs, f3$state_vecs, 1e-3)
})
