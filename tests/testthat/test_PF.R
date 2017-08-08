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
    risk_obj = risk_set,
    F = diag(1, n_vars + 1),
    n_max = 10,
    order = 1,
    n_threads = 1,
    N_fw_n_bw = 10000,
    N_smooth = 1000,
    forward_backward_ESS_threshold = NULL,
    debug = 100) # TODO: remove
  sink()

  expect_equal(
    rep(1, 11), sapply(sapply(result, "[", "weights"), sum),
    check.attributes = FALSE)

  state_est <- t(sapply(result, function(row){
    colSums(t(row$states) * drop(row$weights))
  }))

  # View(cbind(t(result[[2]]$states), result[[2]]$weights))

  matplot(sims$betas, lty = 1, type = "l")
  matplot(state_est, lty = 2, type = "l", add = TRUE)

})

test_that("PF_smooth gives same results", {
  skip_on_cran()

  PF_smooth <- asNamespace("dynamichazard")$PF_smooth

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

  sink("tmp.txt")
  result <- PF_smooth(
    n_fixed_terms_in_state_vec = 0,
    X = t(X_Y$X),
    fixed_terms = t(X_Y$fixed_terms),
    tstart = X_Y$Y[1, ],
    tstop = X_Y$Y[2, ],
    a_0 = sims$betas[1, ],
    Q_0 = diag(1, n_vars + 1),
    Q = diag(1e-2, n_vars + 1),
    risk_obj = risk_set,
    F = diag(1, n_vars + 1),
    n_max = 10,
    order = 1,
    n_threads = 1,
    N_fw_n_bw = 1000,
    N_smooth = 1000,
    forward_backward_ESS_threshold = NULL,
    debug = 100)
  sink()
})
