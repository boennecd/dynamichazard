context("Testing print function")

test_that("Print yields the expected results for object returned by ddhazard", {
  # /w logit
  expect_known_output(
    print(ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 1,
      control = ddhazard_control(eps = 1e-1),
      a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(.1, 2),
      max_T = 20, order = 1)),
    file =  "print_ddhazard", update = FALSE, print = FALSE)

  # /w exponential
  expect_known_output(
    print(ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, model = "exponential",
      control = ddhazard_control(eps = 1e-1),
      a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(.1, 2),
      max_T = 20, order = 1)),
    file =  "print_ddhazard_exp", update = FALSE, print = FALSE)

  # /w UKF
  expect_known_output(
    print(ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, model = "exponential",
      control = ddhazard_control(eps = 1e-1, method = "UKF"),
      a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(.1, 2),
      max_T = 20, order = 1)),
    file =  "print_ddhazard_UKF", update = FALSE, print = FALSE)
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
    control = ddhazard_control(
      fixed_terms_method = "E_step", eps = 1e-2 # decreased to reduce run-time
      ),
    Q_0 = diag(1, 9),
    Q = diag(1e-2, 9),
    max_T = 10, model = "exponential",
    id = sims$res$id, order = 1,
    verbose = F)

  set.seed(999)
  boot_out <- ddhazard_boot(result, R = 19)

  # plot(result, ddhazard_boot = boot_out) # needs bigger R
  expect_known_output(
    print(boot_out, digits = 1), "boot_print", update = FALSE, print = FALSE)

  #####
  # Without fixed effects
  result <- ddhazard(
    survival::Surv(tstart, tstop, event) ~ . - id,
    sims$res,
    by = 1,
    control = ddhazard_control(
      fixed_terms_method = "E_step",
      eps = 1e-2), # decreased to reduce run-time
    Q_0 = diag(1, 11),
    Q = diag(1e-2, 11),
    max_T = 10, model = "exponential",
    id = sims$res$id, order = 1,
    verbose = F)

  set.seed(1992)
  boot_out <- ddhazard_boot(result, R = 19)

  expect_known_output(
    print(boot_out, digits = 1), "boot_print_w_o_fixed", update = FALSE,
    print = FALSE)
})

test_that("Print function for PF objects gives previous results", {
  skip_on_cran()

  .lung <- lung[!is.na(lung$ph.ecog), ]
  set.seed(43588158)
  pf_fit <- suppressWarnings(PF_EM(
    Surv(time, status == 2) ~ ph.ecog + age,
    data = .lung, by = 50, id = 1:nrow(.lung),
    Q_0 = diag(c(1, 1, .1)^2), Q = diag(c(.3, .3, .01)^2),
    max_T = 800,
    control = list(
      N_fw_n_bw = 200,
      N_first = 200,
      N_smooth = 200,
      n_max = 5,
      n_threads = 1)))

  expect_known_output(output <- print(pf_fit), "PF_EM_print")
  expect_equal(output, pf_fit)

  expect_known_output(output <- print(pf_fit$clouds), "PF_clouds_print")
  expect_equal(output, pf_fit$clouds)
})

test_that("Print function for ddhazard_space_errors object gives previous results", {
  fit <- ddhazard(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(eps = 1e-1),
    a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(.1, 2),
    max_T = 20, order = 1)

  errs <- residuals(fit, type = "std_space_error")
  expect_known_output(output <- print(errs), "std_space_error_print")
  expect_equal(output, errs)

  errs <- residuals(fit, type = "space_error")
  expect_known_output(output <- print(errs), "space_error_print")
  expect_equal(output, errs)
})

test_that("print print.ddsurvcurve gives correct result", {
  h1 <- suppressWarnings(ddhazard(
    Surv(stop, event) ~ group, head_neck_cancer, by = 1, max_T = 45,
    Q_0 = diag(2^2, 2), Q = diag(.01^2, 2), control = ddhazard_control(
      method = "GMA", eps = 1e-1, n_max = 1)))
  ddcurve <- ddsurvcurve(h1, new_data = data.frame(
    group = factor(2, levels = 1:2)))

  expect_known_output(ddcurve, file = "print-ddsurvcurve", print = TRUE)
})
