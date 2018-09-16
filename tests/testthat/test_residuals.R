context("Testing residuals functions")

###############
# Simple test that methods calls succeds

cl <- quote(ddhazard(
  formula = survival::Surv(tstart, tstop, event) ~ x1,
  data = logit_sim_200$res,
  id = logit_sim_200$res$id,
  by = 1, max_T = 10,
  a_0 = rep(0, 2), Q_0 = diag(1, 2),
  Q = diag(1e-2, 2)))

cl_logit <- quote(ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ x1 + x2,
    data = logit_sim_200$res,
    by = (by_ <- 1),
    Q_0 = diag(1e-2, 3),
    Q = diag(1000000, 3),
    control = ddhazard_control(est_Q_0 = F, eps = .1), # Just want a fit
    max_T = 10,
    id = logit_sim_200$res$id))

result <- eval(cl)
result_ekf <- result

test_that("Calls to residuals should succed", {
  expect_no_error(residuals(result, "std_space_error"))
  expect_no_error(residuals(result, "space_error"))
  expect_no_error(residuals(result, type = "pearson", data = logit_sim_200))
  expect_no_error(residuals(result, type = "raw", data = logit_sim_200))
})

cl$control <- quote(ddhazard_control(method = "UKF"))
result <- eval(cl)

test_that("residuals functions throws error for some types when method is UKF",{
  expect_error(residuals(result, "std_space_error"))
  expect_error(residuals(result, "space_error"))
  expect_no_error(residuals(result, type = "pearson", data = logit_sim_200))
  expect_no_error(residuals(result, type = "raw", data = logit_sim_200))
})

test_that("Residuals work when data is saved on the fit", {
  for(ty in c("pearson", "raw")){
    cl_cur <- cl_logit
    cl_cur$control <- quote(ddhazard_control(save_data = T))
    fit_saved_data <- eval(cl_cur)

    expect_true(!is.null(fit_saved_data$data))
    fit_saved_data$res <- residuals(fit_saved_data, type = ty)

    cl_cur$control <- quote(ddhazard_control(save_data = F))
    fit_not_saved_data <- eval(cl_cur)

    expect_true(is.null(fit_not_saved_data$data))
    expect_error(residuals(fit_not_saved_data, type = ty))
    fit_not_saved_data$res <- residuals(fit_not_saved_data, type = ty, data = logit_sim_200$res)

    for(i in seq_len(max(length(fit_saved_data$res$residuals),
                         length(fit_not_saved_data$residuals))))
      expect_equal(fit_saved_data$res$residuals[[i]],
                   fit_not_saved_data$res$residuals[[i]])
  }
})

test_that("Residuals work when data is saved on the fit and fixed effects are present", {
  cl_new <- cl_logit
  cl_new$formula <- bquote(update(
    .(cl_new$formula), . ~ . - x1 + ddFixed(x1)))
  cl_new$Q_0 <- diag(1e8, 2)
  cl_new$Q <- diag(.1, 2)

  for(ty in c("pearson", "raw")){
    cl_new$control <- quote(ddhazard_control(save_data = T))
    fit_saved_data <- eval(cl_new)

    expect_true(!is.null(fit_saved_data$data))
    fit_saved_data$res <- residuals(fit_saved_data, type = ty)

    cl_new$control <- quote(ddhazard_control(save_data = F))
    fit_not_saved_data <- eval(cl_new)

    expect_true(is.null(fit_not_saved_data$data))
    expect_error(residuals(fit_not_saved_data, type = ty))
    fit_not_saved_data$res <- residuals(fit_not_saved_data, type = ty, data = logit_sim_200$res)

    for(i in seq_len(max(length(fit_saved_data$res$residuals),
                         length(fit_not_saved_data$residuals))))
      expect_equal(fit_saved_data$res$residuals[[i]],
                   fit_not_saved_data$res$residuals[[i]])
  }
})

########
# Test state space errors

result <-  result_ekf
test_that("State space error match whether standarized or not for logit model", {
  std_res <- residuals(result, "std_space_error")

  expect_true(std_res$standardize)
  expect_s3_class(std_res, "ddhazard_space_errors")

  non_std <- residuals(result, "space_error")
  expect_true(!non_std$standardize)
  expect_s3_class(non_std, "ddhazard_space_errors")

  expect_equal(non_std$Covariances, std_res$Covariances)

  for(i in seq_len(nrow(std_res$residuals)))
    expect_equal(unname(std_res$residuals[i, ]),
                 c(solve(t(chol(non_std$Covariances[, , i]))) %*% non_std$residuals[i, ]))

})

test_that("State space error gives correct dimension with fixed effects", {
  fit <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1, a_0 = 0, Q_0 = as.matrix(1), Q = 1)

  std_res <- residuals(fit, "std_space_error")
  expect_equal(dim(std_res$residuals), c(59, 1))

  suppressWarnings(fit <- ddhazard(
    formula = survival::Surv(start, stop, event) ~ -1 + ddFixed(group),
    data = head_neck_cancer, by = 1))

  expect_warning(std_res <- residuals(fit, "std_space_error"))
  expect_null(std_res)
})


######
# Data sets used in vignettes

# PBC dataset described in Fleming & Harrington (1991)
local({
  pbc <- pbc[!is.na(pbc$protime), ]

  fit <- ddhazard(
    formula = survival::Surv(time, status == 2) ~ log(protime),
    data = pbc, Q_0 = diag(rep(1e8, 2)), by = 100, id = pbc$id,
    Q = diag(rep(1e-2, 2)), max_T = 3600,
    a = rep(0, 2),
    control = ddhazard_control(est_Q_0 = F, eps = .1))

  test_that("Pearson residuals and raw residuals for logistic model are consistent with each other", {
    pearson_res <- residuals(object = fit, type = "pearson", data = pbc)
    raw_res <- residuals(object = fit, type = "raw", data = pbc)

    expect_s3_class(pearson_res, "ddhazard_residual")
    expect_s3_class(raw_res, "ddhazard_residual")

    expect_equal(pearson_res$type, "pearson")
    expect_equal(raw_res$type, "raw")

    for(i in seq_along(pearson_res$residuals))
      expect_equal(pearson_res$residuals[[i]][, "residuals"],
                   raw_res$residuals[[i]][, "residuals"]/
                     sqrt(raw_res$residuals[[i]][, "p_est"] * (1 - raw_res$residuals[[i]][, "p_est"])))
  })

  test_that("Cases in residuals match cases in data", {
    is_case <- pbc$status == 2 & pbc$time <= 3600
    pearson_res <- residuals(object = fit, type = "pearson", data = pbc)

    is_case_residuals <- rep(0, nrow(pbc))
    for(i in seq_along(pearson_res$residuals)){
      y_n_index <- data.frame(pearson_res$residuals[[i]][, c("residuals", "Y", "row_num")])

      is_case_residuals[y_n_index$row_num[y_n_index$Y == 1]] =
        is_case_residuals[y_n_index$row_num[y_n_index$Y == 1]] + 1

      expect_equal(sum(y_n_index$residuals[y_n_index$Y == 1] > 0),
                   sum(y_n_index$Y == 1))
      expect_equal(sum(y_n_index$residuals[y_n_index$Y == 0] < 0),
                   sum(y_n_index$Y == 0))
    }

    expect_equal(is_case_residuals, is_case + 0)
  })
})


##################
# Exponential model

fit <- ddhazard(
  Surv(tstart, tstop, event) ~ x1,
  by = 1, data = exp_sim_200$res,
  id = exp_sim_200$res$id,
  model = "exp_clip_time_w_jump",
  a_0 = numeric(2),
  max_T = 10, Q_0 = diag(1e8, 2), Q = diag(.1, 2),
  control = ddhazard_control(eps = .1))

test_that("Calls to residuals should fail for exponential model and state space error",{
  expect_error(residuals(fit, "std_space_error", regexp = "Functions for with model 'exponential' is not implemented"))
  expect_error(residuals(fit, "space_error", regexp = "Functions for with model 'exponential' is not implemented"))
  expect_error(residuals(fit, "pearson", regexp = "Pearsons residuals is not implemented for model 'exp_clip_time'"))
})

#####
# Test previous result for examples in diagnostics vignette

test_that("Gets previous results with Rossi", {
  skip_on_cran()

  if(interactive()){
    root <- gsub("(.+dynamichazard)(.*)", "\\1", getwd())
    diag_data_path <- paste0(root, "/vignettes/Diagnostics")
  } else{
    diag_data_path <- "."
  }

  load(paste0(diag_data_path, "/Rossi.RData"))

  suppressMessages(
    dd_fit <- dd_fit_first <-  ddhazard(
      Surv(start, stop, event) ~ fin + age + prio + employed.cumsum,
      data = Rossi, id = Rossi$id, by = 1, max_T = 52,
      Q_0 = diag(10000, 5), Q = diag(.1, 5)))

  state <- residuals(dd_fit, type = "std_space_error")
  expect_known_value(state, "Rossi_state.RDS", update = FALSE)

  pearson <- residuals(dd_fit, type = "pearson")
  pearson <- pearson[[1]][1:2]
  expect_known_value(pearson, "Rossi_pearson.RDS", update = FALSE,
                     check.attributes = FALSE)
})

test_that("Get prevoius residuals with whas500", {
  skip_on_cran()

  if(interactive()){
    root <- gsub("(.+dynamichazard)(.*)", "\\1", getwd())
    diag_data_path <- paste0(root, "/vignettes/Diagnostics")
  } else{
    diag_data_path <- "."
  }

  load(paste0(diag_data_path, "/whas500.RData"))

  suppressMessages(
    dd_fit <- ddhazard(
      Surv(lenfol, fstat) ~ gender + age + bmi + hr + cvd,
      data = whas500, by = 100, max_T = 2000,
      Q_0 = diag(10000, 6), Q = diag(.1, 6)))

  state <- residuals(dd_fit, type = "std_space_error")

  expect_known_value(state, "whas_state.RDS", update = FALSE)

  pearson <- residuals(dd_fit, type = "pearson")
  pearson <- pearson[[1]][1:2]

  # save_to_test(pearson, "whas_pearson")
  expect_equal(pearson, read_to_test("whas_pearson"), check.attributes = FALSE)
})
