context("Testing print function")

test_that("Print yields the expected results for object returned by ddhazard", {
  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(1, 2),
    max_T = 20, order = 1)

  # print(paste0(capture.output(print(result, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(result, digits = 4),
              prints_text("Formula:\nsurvival::Surv(stop, event) ~ group\n\nEstimated with EKF in 2 iterations of the EM algorithm\n\nEstimated time-varying effects and point-wise standard deviation:\n    (Intercept)    sd   group1    sd \nt0       -1.802 0.6527 0.04366 0.7072\nt1       -2.337 0.2933 0.08239 0.3937\nt2       -2.803 0.3769 0.04722 0.4965\nt3       -2.786 0.4343 0.47588 0.5711\nt4       -2.434 0.4577 0.03763 0.5565\nt5       -1.785 0.4261 0.80028 0.5652\nt6       -1.973 0.3738 0.66434 0.4719\nt7       -2.450 0.4073 0.19786 0.5188\nt8       -2.878 0.4576 0.70849 0.6143\nt9       -3.016 0.5184 0.74029 0.6516\nt10      -3.094 0.5524 0.72503 0.6899\nt11      -3.389 0.5742 0.40016 0.7260\nt12      -3.600 0.6127 0.44614 0.8121\nt13      -3.734 0.6424 0.99283 0.8730\nt14      -3.113 0.6846 2.14297 0.8807\nt15      -3.279 0.6211 1.62011 0.7171\nt16      -3.576 0.6459 1.14035 0.7634\nt17      -3.895 0.6823 1.04613 0.8475\nt18      -3.689 0.7507 1.36408 0.9465\nt19      -3.951 0.7803 1.55400 1.0122\nt20      -3.613 0.9415 2.34599 1.2115",
                          fixed = T))

  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ ddFixed_intercept() + group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = F, save_risk_set = F,
                   fixed_terms_method = "E_step"),
    a_0 = rep(0, 1), Q_0 = diag(1, 1), Q = 1,
    max_T = 20, order = 1)

  # print(paste0(capture.output(print(result, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(result, digits = 4),
              prints_text("Formula:\nsurvival::Surv(stop, event) ~ ddFixed_intercept() + group\n\nEstimated with EKF in 2 iterations of the EM algorithm\n\nEstimated time-varying effects and point-wise standard deviation:\n      group1    sd \nt0  -0.30363 0.7490\nt1  -0.40003 0.5721\nt2   0.14752 0.6353\nt3   0.89512 0.5788\nt4   0.61924 0.4402\nt5   1.71482 0.5019\nt6   1.52557 0.4088\nt7   0.66129 0.4474\nt8   0.84395 0.5513\nt9   0.69078 0.5495\nt10  0.53181 0.5825\nt11 -0.01786 0.6187\nt12 -0.14110 0.7235\nt13  0.31465 0.7936\nt14  1.52097 0.7866\nt15  1.02346 0.5449\nt16  0.38831 0.6161\nt17  0.16822 0.7148\nt18  0.49660 0.8016\nt19  0.65475 0.8562\nt20  1.53209 1.0514\n\nFixed effects are estimated in the E_step. The estimates are:\n(Intercept) \n     -2.926 ",
                          fixed = T))

  suppressWarnings(result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ ddFixed_intercept() + ddFixed(as.numeric(group)),
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = F, save_risk_set = F,
                   fixed_terms_method = "E_step"),
    max_T = 20, order = 1))

  # print(paste0(capture.output(print(result, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(result, digits = 4),
              prints_text("Formula:\nsurvival::Surv(stop, event) ~ ddFixed_intercept() + ddFixed(as.numeric(group))\n\nEstimated with EKF in 3 iterations of the EM algorithm\n\nFixed effects are estimated in the E_step. The estimates are:\n      (Intercept) as.numeric(group) \n          -3.4628            0.5524 ",
                          fixed = T))
})

test_that("print.ddhazard_boot gives the expected output", {
  skip_on_cran()

  #####
  # With fixed effects

  sims <- exp_sim_200
  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ ddFixed(x1) + ddFixed(x2) + . - id - x1 - x2,
    sims$res,
    by = 1,
    control = list(fixed_terms_method = "E_step",
                   eps = 1e-2), # decreased to reduce run-time
    Q_0 = diag(1, 9),
    Q = diag(1e-2, 9),
    max_T = 10, model = "exp_clip_time_w_jump",
    id = sims$res$id, order = 1,
    verbose = F)

  set.seed(999)
  boot_out <- ddhazard_boot(result, R = 19)

  # plot(result, ddhazard_boot = boot_out) # needs bigger R
  # save_to_test(
  #   capture.output(print(boot_out, digits = 1)),
  #   file_name = "boot_print.RDS")
  expect_equal(capture.output(print(boot_out, digits = 1)), read_to_test("boot_print.RDS"))

  #####
  # Without fixed effects
  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ . - id,
    sims$res,
    by = 1,
    control = list(fixed_terms_method = "E_step",
                   eps = 1e-2), # decreased to reduce run-time
    Q_0 = diag(1, 11),
    Q = diag(1e-2, 11),
    max_T = 10, model = "exp_clip_time_w_jump",
    id = sims$res$id, order = 1,
    verbose = F)

  set.seed(1992)
  boot_out <- ddhazard_boot(result, R = 19)

  # save_to_test(
  #   capture.output(print(boot_out, digits = 1)),
  #   file_name = "boot_print_two.RDS")
  expect_equal(capture.output(print(boot_out, digits = 1)), read_to_test("boot_print_two.RDS"))
})
