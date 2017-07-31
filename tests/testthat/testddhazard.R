context("Testing ddhazard function")

# Test on data set that is one of Farhmiers papers
result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  control = list(est_Q_0 = F,
                 save_data = F, save_risk_set = F),
  a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
  max_T = 45,
  id = head_neck_cancer$id, order = 1)

# plot(result)

test_that("Testing names of output from ddhazard on head and neck cancer dataset", {
  expect_equal(colnames(result$state_vecs), c("(Intercept)", "group1"))
  expect_equal(unlist(dimnames(result$state_vars)), unlist(list(c("(Intercept)", "group1"), c("(Intercept)", "group1"), NULL)))
  expect_equal(unlist(dimnames(result$Q)), rep(c("(Intercept)", "group1"), 2))
  expect_equal(unlist(dimnames(result$Q_0)), rep(c("(Intercept)", "group1"), 2))
})

result <- result[c("state_vecs", "state_vars","lag_one_cov", "model")]
# save_to_test(result, "ddhazard_head_neck")

test_that("get previous results with head_neck", {
  expect_equal(result, read_to_test("ddhazard_head_neck"))
})



test_that("Invalid penalty terms throw error", {
  expect_error(
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, # Use by month intervals
      control = list(denom_term = 0)),
    regexp = "Method not implemented with penalty term \\(control\\$denom_term\\) equal to 0")

  expect_error(
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, # Use by month intervals
      control = list(denom_term = -1)),
    regexp = "Method not implemented with penalty term \\(control\\$denom_term\\) equal to -1")
})



test_that("Changing convergence criteria change output",{
  arg_list <- list(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, # Use by month intervals
    id = head_neck_cancer$id,
    Q_0 = diag(1e5, 2), Q = diag(.1, 2),
    control = list(criteria = "delta_coef", eps = .002))

  suppressMessages(res1 <- do.call(ddhazard, arg_list))
  arg_list$control$criteria <- "delta_likeli"
  suppressMessages(res2 <- do.call(ddhazard, arg_list))

  expect_true(res1$n_iter != res2$n_iter)
})

suppressMessages(
  result_exp <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 3, Q_0 = diag(10000, 2),
    Q = diag(1e-3, 2),
    control = list(est_Q_0 = F, n_max = 10^3, eps = 10^-4,
                   denom_term = 6e-2, save_data = F, save_risk_set = F),
    max_T = 30,
    id = head_neck_cancer$id, order = 1,
    verbose = F,
    model = "exp_clip_time_w_jump"))

test_that("Result of exponential model on head_neck_data match previous results", {
  # plot(result_exp)
  result_exp$control <- NULL
  result_exp$call <- NULL
  # save_to_test(result_exp, "head_neck_exp")

  expect_equal(result_exp, read_to_test("head_neck_exp"))
})

test_that("exponential model and logit moels hazzard functions differs", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F,
                   save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  expect_true(result_exp$model != result$model)
  expect_true(toString(body(result_exp$hazard_func)) !=
                toString(body(result$hazard_func)))
  expect_true(toString(body(result_exp$hazard_first_deriv)) !=
                toString(body(result$hazard_first_deriv)))
})

test_that("Testing names of output from ddhazard on head and neck cancer dataset", {
  expect_equal(colnames(result_exp$state_vecs), c("(Intercept)", "group1"))
  expect_equal(unlist(dimnames(result_exp$state_vars)), unlist(list(c("(Intercept)", "group1"), c("(Intercept)", "group1"), NULL)))
  expect_equal(unlist(dimnames(result_exp$Q)), rep(c("(Intercept)", "group1"), 2))
  expect_equal(unlist(dimnames(result_exp$Q_0)), rep(c("(Intercept)", "group1"), 2))
})

test_that("You can ommit the first entry when covariates are not time-varying",{
  r1 <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, # Use by month intervals
    control = list(est_Q_0 = F, n_max = 10^4, eps = 10^-4,
                   save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(1, 2), # Initial value
    max_T = 45,
    id = head_neck_cancer$id, order = 1
  )

  r2 <- ddhazard(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, # Use by month intervals
    control = list(est_Q_0 = F, n_max = 10^4, eps = 10^-4,
                   save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(1, 2), # Initial value
    max_T = 45,
    id = head_neck_cancer$id, order = 1
  )

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_equal(r1$state_vars, r2$state_vars)
})


# Change by argument
test_that("Chaning by argument gives previous results for the exp_clp_w_jump method", {
  suppressMessages(result_exp <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 5, Q_0 = diag(1, 2),
    Q = diag(1e-1, 2),
    control = list(est_Q_0 = F,
                   save_data = F, save_risk_set = F),
    max_T = 30,
    id = head_neck_cancer$id, order = 1,
    verbose = F,
    model = "exp_clip_time_w_jump"
  ))

  # plot(result_exp)
  result_exp$control <- NULL
  result_exp$call <- NULL
  # save_to_test(result_exp, "ddhazard_changed_by")

  expect_equal(result_exp, read_to_test("ddhazard_changed_by"))
})

test_that("Unmacthed control variable throw error",
          expect_error({
            result = ddhazard(
              formula = survival::Surv(start, stop, event) ~ group,
              data = head_neck_cancer,
              by = 1, # Use by month intervals
              a_0 = rep(0, 2), Q_0 = diag(1, 2), # Initial value
              max_T = 45,
              id = head_neck_cancer$id, order = 1,
              control = list(None_existing_parem = 1)
            )}, regexp = "These control parameters are not recognized"))

########
# Test on simulated data

test_that("Result of exponential model with only binary or right clipped time yield previous results", {
  args <- list(
    formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    data = exp_sim_200$res,
    by = 1,
    Q_0 = diag(10000, 11),
    Q = diag(1e-2, 11),
    control = list(
      est_Q_0 = F, eps = 10^-2,
      save_data = F, save_risk_set = F,
      method = "EKF"),
    max_T = 10,
    id = exp_sim_200$res$id, order = 1,
    model = "exp_bin")

  result_exp <- do.call(ddhazard, args)

  # matplot(exp_sim_200$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, lty = 2, type = "l", add = T)
  result_exp <- result_exp[c("state_vars", "state_vecs", "Q")]
  # save_to_test(result_exp, "ddhazard_exp_bin")

  expect_equal(result_exp, read_to_test("ddhazard_exp_bin"), tolerance = 1e-5)

  args$model <- "exp_clip_time"
  result_exp <- do.call(ddhazard, args)

  # matplot(exp_sim_200$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, lty = 2, type = "l", add = T)

  result_exp <- result_exp[c("state_vecs", "state_vars", "Q")]
  # save_to_test(result_exp, file_name = "ddhazard_exp_clip", tolerance = 1e-4)
  expect_equal(result_exp, read_to_test("ddhazard_exp_clip"), tolerance = 1e-04)

  args$model <- "exp_clip_time_w_jump"
  result_exp <- do.call(ddhazard, args)

  # matplot(exp_sim_200$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, lty = 2, type = "l", add = T)

  result_exp <- result_exp[c("state_vecs", "state_vars", "Q")]
  # save_to_test(result_exp, file_name = "ddhazard_exp_clip_w_jump", 1e-4)
  expect_equal(result_exp, read_to_test("ddhazard_exp_clip_w_jump"), tolerance = 1e-04)
})


test_that("Permutating data does not makes a big difference", {
  args <- list(
    Surv(stop, event) ~ group, head_neck_cancer,
    by = 1, max_T = 40,
    Q_0 = diag(rep(10000, 2)), Q = diag(rep(0.1, 2)),
    control = list(n_threads = 1))

  r1 <- do.call(ddhazard, args)

  args$control <- c(args$control, list(permu = T))

  r2 <- do.call(ddhazard, args)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(is.character(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))

  #####
  # With fixed effects
  args <- list(
    Surv(tstart, tstop, death == 2) ~ age + ddFixed(edema) +
      log(albumin) + log(protime) + log(bili), pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    Q_0 = diag(rep(10000, 5)), Q = diag(rep(0.001, 5)),
    control = list(n_threads = 1))

  r1 <- do.call(ddhazard, args)
  args$control <- c(args$control, list(permu = T))
  r2 <- do.call(ddhazard, args)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(is.character(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))


  #####
  # With weigths
  w <- sample(1:3, nrow(pbc2), replace = T)

  args <- list(
    Surv(tstart, tstop, death == 2) ~ age + edema +
      log(albumin) + log(protime) + log(bili), pbc2,
    id = pbc2$id, by = 100, max_T = 3000,
    weights = w,
    Q_0 = diag(rep(10000, 6)), Q = diag(rep(0.001, 6)),
    control = list(n_threads = 1, LR = .8))

  r1 <- do.call(ddhazard, args)
  args$control <- c(args$control, list(permu = T))
  r2 <- do.call(ddhazard, args)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(is.character(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))

})


