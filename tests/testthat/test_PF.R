context("Running test_PF")

test_that("dmvnrm_log_test gives correct likelihood", {
  skip_if_not_installed("mvtnorm")

  dmvnrm_log_test <- asNamespace("dynamichazard")$dmvnrm_log_test

  for(i in 1:10){
    n <- 5
    mu <- rnorm(n)
    A <- matrix(runif(n^2)*2-1, ncol=n)
    sigma <- t(A) %*% A
    sigma_chol_inv <- solve(chol(sigma))

    for(i in 1:10){
      x <- rnorm(n)
      expect_equal(
        mvtnorm::dmvnorm(x, mu, sigma, log = TRUE),
        dmvnrm_log_test(x, mu, sigma_chol_inv))
    }
  }
})

test_that("FW_filter gives gives same results", {
  skip_on_cran()

  FW_filter <- asNamespace("dynamichazard")$FW_filter

  n_vars <- 2
  set.seed(78095324)
  sims <- test_sim_func_logit(
    n_series = 100, n_vars = n_vars, t_0 = 0, t_max = 10,
    x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
    intercept_start = -3, sds = c(.1, rep(.5, n_vars)))

  # sum(sims$res$event)
  # matplot(sims$beta, type = "l", lty = 1)

  X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)

  risk_set <-
    get_risk_obj(Y = X_Y$Y, by = 1, max_T = 10,
                 id = sims$res$id, is_for_discrete_model = TRUE)

  old_seed <- .Random.seed
  sink("tmp.txt")
  result <- FW_filter(
    n_fixed_terms_in_state_vec = 0,
    X = t(X_Y$X),
    fixed_terms = t(X_Y$fixed_terms),
    tstart = X_Y$Y[1, ],
    tstop = X_Y$Y[2, ],
    a_0 = sims$betas[1, ],
    Q_0 = diag(1, n_vars + 1),
    Q = diag(1, n_vars + 1),
    Q_tilde = diag(1, n_vars + 1),
    risk_obj = risk_set,
    F = diag(1, n_vars + 1),
    n_max = 10,
    order = 1,
    n_threads = 1,
    N_fw_n_bw = 1000,
    N_smooth = 1000,
    forward_backward_ESS_threshold = NULL,
    debug = 100) # TODO: remove
  sink()

  expect_false(all(old_seed == .Random.seed))

  expect_equal(
    rep(1, 11), sapply(sapply(result, "[", "weights"), sum),
    check.attributes = FALSE)

  state_est <- t(sapply(result, function(row){
    colSums(t(row$states) * drop(row$weights))
  }))

  # View(cbind(t(result[[2]]$states), result[[2]]$weights))

  matplot(sims$betas, lty = 1, type = "l")
  matplot(state_est, lty = 2, type = "l", add = TRUE)

  exect_true(FALSE)
})






test_that("BW_filter gives gives same results", {
  skip_on_cran()

  BW_filter <- asNamespace("dynamichazard")$BW_filter

  n_vars <- 2
  set.seed(78095324)
  sims <- test_sim_func_logit(
    n_series = 100, n_vars = n_vars, t_0 = 0, t_max = 10,
    x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
    intercept_start = -3, sds = c(.1, rep(.5, n_vars)))

  # sum(sims$res$event)
  # matplot(sims$beta, type = "l", lty = 1)

  X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)

  risk_set <-
    get_risk_obj(Y = X_Y$Y, by = 1, max_T = 10,
                 id = sims$res$id, is_for_discrete_model = TRUE)

  old_seed <- .Random.seed
  sink("tmp.txt")
  result <- BW_filter(
    n_fixed_terms_in_state_vec = 0,
    X = t(X_Y$X),
    fixed_terms = t(X_Y$fixed_terms),
    tstart = X_Y$Y[1, ],
    tstop = X_Y$Y[2, ],
    a_0 = sims$betas[1, ],
    Q_0 = diag(1, n_vars + 1),
    Q = diag(1, n_vars + 1),
    Q_tilde = diag(1, n_vars + 1),
    risk_obj = risk_set,
    F = diag(1, n_vars + 1),
    n_max = 10,
    order = 1,
    n_threads = 1,
    N_fw_n_bw = 1000,
    N_smooth = 1000,
    forward_backward_ESS_threshold = NULL,
    debug = 100) # TODO: remove
  sink()

  expect_false(all(old_seed == .Random.seed))

  expect_equal(
    rep(1, 11), sapply(sapply(result, "[", "weights"), sum),
    check.attributes = FALSE)

  state_est <- t(sapply(result, function(row){
    colSums(t(row$states) * drop(row$weights))
  }))

  state_est <- state_est[11:2, ]

  # View(cbind(t(result[[2]]$states), result[[2]]$weights))

  matplot(sims$betas, lty = 1, type = "l", ylim = range(sims$betas, state_est))
  matplot(2:11, state_est, lty = 2, type = "l", add = TRUE)

  exect_true(FALSE)
})





test_that("PF_smooth gives same results", {
  skip_on_cran()

  PF_smooth <- asNamespace("dynamichazard")$PF_smooth

  n_vars <- 2
  set.seed(78095324)
  sims <- test_sim_func_logit(
    n_series = 250, n_vars = n_vars, t_0 = 0, t_max = 10,
    x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
    intercept_start = -3, sds = c(.1, rep(.5, n_vars)))

  # sum(sims$res$event)
  # matplot(sims$beta, type = "l", lty = 1)

  X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)

  risk_set <-
    get_risk_obj(Y = X_Y$Y, by = 1, max_T = 10,
                 id = sims$res$id, is_for_discrete_model = TRUE)

  old_seed <- .Random.seed
  sink("tmp.txt")
  result <- PF_smooth(
    n_fixed_terms_in_state_vec = 0,
    X = t(X_Y$X),
    fixed_terms = t(X_Y$fixed_terms),
    tstart = X_Y$Y[1, ],
    tstop = X_Y$Y[2, ],
    a_0 = sims$betas[1, ],
    Q_0 = diag(1, n_vars + 1),
    Q = diag(1, n_vars + 1),
    Q_tilde = diag(1, n_vars + 1),
    risk_obj = risk_set,
    F = diag(1, n_vars + 1),
    n_max = 10,
    order = 1,
    n_threads = 1,
    N_fw_n_bw = 1000,
    N_smooth = 10000,
    forward_backward_ESS_threshold = NULL,
    debug = 2)
  sink()

  expect_false(all(old_seed == .Random.seed))

  matplot(0:10, sims$betas, lty = 1, type = "l", ylim = c(-5, 5), xlim = c(0, 11))
  for(i in 1:3){
    state_est <- t(sapply(result[[i]], function(row){
      colSums(t(row$states) * drop(row$weights))
    }))

    idx <- switch(i, "1" = 0:10, "2" = 1:11, "3" = 1:10)
    matplot(idx, state_est, lty = i + 1, type = "l", add = TRUE)
    matplot(idx, state_est, lty = i + 1, type = "p", add = TRUE, pch = 16 + i)
  }

  sapply(lapply(result$smoothed_clouds, "[[", "parent_idx"),
         function(n) sort(xtabs(~ n), decreasing = TRUE)[1:10])
  sapply(lapply(result$smoothed_clouds, "[[", "parent_idx"),
         function(n) sort(xtabs(~ n), decreasing = TRUE)[1:10])

  result <- result$smoothed_clouds

  expect_equal(
    rep(1, 10), sapply(sapply(result, "[", "weights"), sum),
    check.attributes = FALSE)

  exect_true(FALSE)
})
