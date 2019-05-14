context("testing cloglog link function")

set.seed(57826639)
dat <- test_sim_func_discrete(
  n_series = 400, n_vars = 2, t_0 = 0, t_max = 30,
  x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(2),
  is_fixed = 2,
  intercept_start = -4, sds = c(.1, rep(.5, 2)), linkfunc = "cloglog")

test_that("cloglog function gives previous result with simulated data using ddhazard", {
  fit <- ddhazard(
    Surv(tstart, tstop, event) ~ ddFixed(x1) + x2, data = dat$res, by = 1L,
    id = dat$res$id, Q = diag(1, 2), Q_0 = diag(1, 2), model = "cloglog",
    control = ddhazard_control(method = "EKF"))

  # matplot(fit$state_vecs, type = "l", lty = 2)
  # matlines(dat$betas[, c(1, 3)], type = "l", lty = 1)

  do_check <- c("state_vecs", "state_vecs", "Q", "a_0", "fixed_effects")
  expect_known_value(
    fit[names(fit) %in% do_check], file = "cloglog-sim-ekf-1-it.RDS", update = FALSE,
    check.attributes = FALSE)

  fit <- ddhazard(
    Surv(tstart, tstop, event) ~ ddFixed(x1) + x2, data = dat$res, by = 1L,
    id = dat$res$id, Q = diag(1, 2), Q_0 = diag(1, 2), model = "cloglog",
    control = ddhazard_control(method = "EKF", NR_eps = 1e-4))
  # matplot(fit$state_vecs, type = "l", lty = 2)
  # matlines(dat$betas[, c(1, 3)], type = "l", lty = 1)

  expect_known_value(
    fit[names(fit) %in% do_check], file = "cloglog-sim-ekf-more-it.RDS", update = FALSE,
    check.attributes = FALSE)

  fit <- ddhazard(
    Surv(tstart, tstop, event) ~ ddFixed(x1) + x2, data = dat$res, by = 1L,
    id = dat$res$id, Q = diag(1, 2), Q_0 = diag(1, 2), model = "cloglog",
    control = ddhazard_control(method = "GMA"))
  # matplot(fit$state_vecs, type = "l", lty = 2)
  # matlines(dat$betas[, c(1, 3)], type = "l", lty = 1)
  expect_known_value(
    fit[names(fit) %in% do_check], file = "cloglog-sim-gma.RDS", update = FALSE,
    check.attributes = FALSE)

  preds <- predict(
    fit, new_data = data.frame(
      tstart = c(0, 0), tstop = c(30, 30),
      x1 = c(-.5, .5), x2 = c(-.5, .5)),
    type = "response", tstart = "tstart", tstop = "tstop")
  expect_known_value(
    preds, file = "cloglog-preds.RDS", update = FALSE,
    check.attributes = FALSE)
})

test_that("cloglog function gives previous result with simulated data using PF_EM and PF_forward_filter", {
  skip_on_cran()
  fi <-file.path("local_tests", "cloglog-PF-RW.RDS")
  skip_if(!file.exists(file.path("previous_results", fi)))

  set.seed(1)
  fit <- suppressWarnings( # warns that method did not converge
    PF_EM(
      Surv(tstart, tstop, event) ~ ddFixed(x1) + x2, data = dat$res, by = 1L,
      id = dat$res$id, Q = diag(1, 2), Q_0 = diag(1, 2), model = "cloglog",
      type = "RW",
      control = PF_control(
        method = "AUX_normal_approx_w_cloud_mean", N_fw_n_bw = 200,
        N_smooth = 500, n_threads = parallel::detectCores(logical = FALSE),
        N_first = 1000L, n_max = 1L, eps = 1e-4)))

  expect_known_value(
    fit[c("EM_ests", "effective_sample_size", "log_ligkes")],
    file = fi, check.attributes = FALSE)

  set.seed(1)
  fw_res <- PF_forward_filter(
    dat$res, N_fw = 100L, N_first = 100L, by = 1L, type = "RW",
    id = dat$res$id, a_0 = c(-4.08959, x2 = 1.8239682), Q_0 = diag(1, 2),
    Q = matrix(c(0.0136552, -0.02407266, -0.02407266, 0.08846885), 2L),
    fixed_effects = 0.75606031, formula =
      Surv(tstart, tstop, event) ~ ddFixed(x1) + x2, control = PF_control(
        N_first = 1L, N_fw_n_bw = 1L, N_smooth = 1L, method =
          "AUX_normal_approx_w_cloud_mean",
        n_threads = parallel::detectCores(logical = FALSE)))

  expect_equal(as.numeric(logLik(fw_res)), -1002.0181117173082157)
})
