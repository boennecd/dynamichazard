context("Testing examples in man files")

test_that("ddhazard help page examples gives the same results", {
  pbc <- pbc_org
  fit <- ddhazard(
    Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
    Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
    control = list(method = "GMA"))

  # test for plot
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, cov_index = 2))
  # end test for plot

  # test for predict
  pred <- predict(fit, type = "response", new_data =
                    data.frame(time = 0, status = 2, bili = 3))
  expect_known_value(pred, "ddhazard_help_page_predict_response.RDS")
  pred <- predict(fit, type = "term", new_data =
                    data.frame(time = 0, status = 2, bili = 3))
  expect_known_value(pred, "ddhazard_help_page_predict_term.RDS")
  # end test for predict

  # test for loglike
  ll <- logLik(fit)
  expect_known_value(ll, "ddhazard_help_page_loglike.RDS")
  # end test for loglike

  fit <- fit[c("state_vecs", "state_vars", "lag_one_cov")]
  expect_known_value(fit, "ddhazard_help_page_order_one.RDS")

  # second order model
  fit <- ddhazard(
   Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
   Q_0 = diag(1, 4), Q = diag(1e-4, 2), by = 50,
   control = list(method = "GMA"),
   order = 2)
  fit <- fit[c("state_vecs", "state_vars", "lag_one_cov")]
  expect_known_value(fit, "ddhazard_help_page_order_two.RDS")
})

test_that("residuals.ddhazard help page examples gives the same results", {
  skip_on_cran()

  pbc <- pbc_org
  fit <- ddhazard(
    Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
    Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
    control = list(method = "GMA"))

  resids <- residuals(fit, type = "pearson")$residuals
  expect_known_value(resids[1:2], "ddhazard_help_page_residuals.RDS")
})

test_that("hatvalues.ddhazard help page examples gives the same results", {
  pbc <- pbc_org
  fit <- ddhazard(
    Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3000,
    Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 100,
    control = list(method = "GMA"))
  hvs <- hatvalues(fit)

  expect_known_value(hvs[1:2], "hatvalues_help_page.RDS")
})

test_that("ddhazard_boot help page examples gives the same results", {
  skip_on_cran()

  pbc <- pbc_org
  set.seed(56219373)
  fit <- ddhazard(
    Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3000,
    Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 100,
    control = list(method = "GMA"))
  bt <- ddhazard_boot(fit, R = 999)

  expect_known_value(
    bt$t[1:10, ], "ddhazard_boot_help_page.RDS", tolerance = 1e-4)
  expect_no_error(suppressMessages(plot(fit, ddhazard_boot = bt, level = .9)))
})

test_that("get_survival_case_weights_and_data help page examples gives the same results", {
  dat <- data.frame(
    id     = c(   1,    1, 2,     2),
    tstart = c(   0,    4, 0,     2),
    tstop  = c(   4,    6, 2,     6),
    event  = c(   0,    1, 0,     0),
    x1     = c(1.09, 1.29, 0, -1.16))

  res <- get_survival_case_weights_and_data(
    Surv(tstart, tstop, event) ~ x1, dat, by = 1, id = dat$id)$X
  expect_known_value(res, "get_survival_case_weights_and_data_help_page_1.RDS")
  res <- get_survival_case_weights_and_data(
    Surv(tstart, tstop, event) ~ x1, dat, by = 1, id = dat$id,
    use_weights = FALSE)$X
  expect_known_value(res, "get_survival_case_weights_and_data_help_page_2.RDS")
})

test_that("get_risk_obj help page examples gives the same results", {
  dat <- data.frame(
    id     = c(1, 1, 2, 2),
    tstart = c(0, 4, 0, 2),
    tstop  = c(4, 6, 2, 4),
    event  = c(0, 1, 0, 0))

  out <- with(dat, get_risk_obj(Surv(tstart, tstop, event), by = 1, max_T = 6, id = id))
  expect_known_value(out, "get_risk_obj_help_page.RDS")
})

test_that("static_glm help page examples gives the same results", {
  pbc <- pbc_org

  fit <- static_glm(
   Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
   by = 50)
  expect_known_value(fit$coefficients, "static_glm_help_page.RDS")
})

# Function to compute cloud means
get_means <- function(result){
  means <- lapply(
    result[c("forward_clouds", "backward_clouds", "smoothed_clouds")],
    function(clouds)
      do.call(rbind, sapply(clouds, function(row){
        colSums(t(row$states) * drop(row$weights))
      }, simplify = FALSE)))
}

test_that("PF_EM help page example runs and gives previous computed results", {
  skip_on_cran()

  .lung <- lung[!is.na(lung$ph.ecog), ]
  set.seed(43588155)
  pf_fit <- PF_EM(
    Surv(time, status == 2) ~ ph.ecog + age,
    data = .lung, by = 50, id = 1:nrow(.lung),
    Q_0 = diag(1, 3), Q = diag(1, 3),
    max_T = 800,
    control = PF_control(
      N_fw_n_bw = 500,
      N_first = 2500,
      N_smooth = 2500,
      n_max = 50,
      n_threads = max(parallel::detectCores(), 2)))

  .file <- "local_tests/survival_lung_example"
  # save_to_test(pf_fit[names(pf_fit) != "call"], .file)
  test_if_file_exists(
    .file,
    expect_equal(pf_fit[names(pf_fit) != "call"], read_to_test(.file),
                 tolerance = 1.49e-08))

  .file <- "survival_lung_example_cloud_means"
  # save_to_test(get_means(pf_fit$clouds), file_name = .file)
  expect_equal(get_means(pf_fit$clouds), read_to_test(.file), tolerance = 1.49e-08)

  expect_no_error(plot(pf_fit, cov_index = 1))
  expect_no_error(plot(pf_fit, cov_index = 2))
  expect_no_error(plot(pf_fit, cov_index = 3))
  expect_no_error(plot(pf_fit$log_likes))
})

test_that("Second example on PF help page gives the same result", {
  skip_on_cran()
  skip_if(!dir.exists("previous_results/local_tests"))

  # prepare data
  pbc <- survival::pbc
  pbcseq <- survival::pbcseq
  temp <- subset(pbc, id <= 312, select=c(id, sex, time, status, edema, age))
  pbc2 <- tmerge(temp, temp, id=id, death = event(time, status))
  pbc2 <- tmerge(pbc2, pbcseq, id=id, albumin = tdc(day, albumin),
                 protime = tdc(day, protime), bili = tdc(day, bili))
  pbc2 <- pbc2[, c("id", "tstart", "tstop", "death", "sex", "edema",
                   "age", "albumin", "protime", "bili")]
  pbc2 <- within(pbc2, {
    log_albumin <- log(albumin)
    log_protime <- log(protime)
    log_bili <- log(bili)
  })

  # standardize
  for(c. in c("age", "log_albumin", "log_protime", "log_bili"))
    pbc2[[c.]] <- drop(scale(pbc2[[c.]]))

  # fit model with extended Kalman filter
  ddfit <- ddhazard(
    Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() + ddFixed(age) +
      ddFixed(edema) + ddFixed(log_albumin) + ddFixed(log_protime) + log_bili,
    pbc2, Q_0 = 100, Q = 1e-2, by = 100, id = pbc2$id,
    model = "exponential", max_T = 3600,
    control = list(eps = 1e-5, NR_eps = 1e-4, n_max = 1e4))

  expect_known_value(ddfit[names(ddfit) != "call"],
                     "local_tests/pf_man_2nd_ddfit.RDS")

  # fit model with particle filter
  set.seed(88235076)
  ppfit <- suppressWarnings(PF_EM(
    Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() + ddFixed(age) +
      ddFixed(edema) + ddFixed(log_albumin) + ddFixed(log_protime) + log_bili,
    pbc2, Q_0 = 100, Q = ddfit$Q * 100, # use estimate from before
    by = 100, id = pbc2$id,
    model = "exponential", max_T = 3600,
    control = PF_control(
      N_fw_n_bw = 250, N_smooth = 500, N_first = 1000, eps = 1e-3,
      method = "AUX_normal_approx_w_cloud_mean",
      n_max = 25, # just take a few iterations as an example
      n_threads = parallel::detectCores() - 2L)))

  expect_known_value(ppfit[!names(ppfit) %in%
                             c("clouds", "call", "summary_stats")],
                     "local_tests/pf_man_2nd_ppfit.RDS")
})
