context("Testing summary function")

test_that("summary.ddhazard yields previous results and prints match previous values", {
  test_func <- function(fit, file, update, ...){
    .summary <- do.call(summary, c(list(object = fit), list(...)))

    expect_known_value(.summary, paste0(file, "_object.RDS"), update = update)
    expect_known_output(
      print(.summary), paste0(file, "_print"), update = update, print = FALSE)

    invisible(list(summary = .summary, fit = fit))
  }

  # /w logit
  fit <- ddhazard(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = ddhazard_control(eps = 1e-1),
    a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(.1, 2),
    max_T = 20, order = 1)
  f1 <- test_func(fit, file =  "summary_ddhazard", update = FALSE)

  # /w logit & all coefs and time
  f2 <- test_func(
    fit, file =  "summary_ddhazard_all_coef", update = FALSE, max_print = Inf)
  expect_true(!isTRUE(all.equal(f1$summary, f2$summary)))

  # /w logit & only one coef
  f3 <- test_func(
    fit, file =  "summary_ddhazard_one_coef", update = FALSE, var_indices = 1)
  expect_true(!isTRUE(all.equal(f2$summary, f3$summary)))

  # /w logit & w/ changed digits
  old <- getOption("digits")
  options(digits = 1)
  f4 <- test_func(
    fit, file =  "summary_ddhazard_one_digit", update = FALSE)
  options(digits = old)
  expect_equal(f1$summary, f4$summary)

  # /w logit and fixed
  test_func(
    ddhazard(
      formula = survival::Surv(stop, event) ~ ddFixed(group),
      data = head_neck_cancer,
      by = 1,
      control = ddhazard_control(eps = 1e-1),
      a_0 = 0, Q_0 = 1, Q = .1,
      max_T = 20, order = 1),
    file =  "summary_ddhazard_fixed", update = FALSE)

  # /w exponential
  test_func(
    ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, model = "exponential",
      control = ddhazard_control(eps = 1e-1),
      a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(.1, 2),
      max_T = 20, order = 1),
    file =  "summary_ddhazard_exp", update = FALSE)

  # /w UKF
  test_func(
    ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, model = "exponential",
      control = ddhazard_control(eps = 1e-1, method = "UKF"),
      a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(.1, 2),
      max_T = 20, order = 1),
    file =  "summary_ddhazard_UKF", update = FALSE)
})
