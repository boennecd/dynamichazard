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


test_that("PF_smooth gives same results", {
  skip_on_cran()

  PF_smooth <- asNamespace("dynamichazard")$PF_smooth

  n_vars <- 2
  set.seed(78095324)
  sims <- test_sim_func_logit(
    n_series = 1e3, n_vars = n_vars, t_0 = 0, t_max = 10,
    x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
    intercept_start = -3, sds = c(.1, rep(.5, n_vars)))

  X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)
  risk_set <-
    get_risk_obj(Y = X_Y$Y, by = 1, max_T = 10,
                 id = sims$res$id, is_for_discrete_model = TRUE)

  # sum(sims$res$event)
  # matplot(sims$beta, type = "l", lty = 1)
  # xtabs(~ sims$res$tstop[sims$res$event == 1])

  Q <- diag(sqrt(.33), 3)
  Q_0 <- diag(.5, 3)
  a_0 <- sims$betas[1, ]

  args <- list(
    n_fixed_terms_in_state_vec = 0,
    X = t(X_Y$X),
    fixed_terms = t(X_Y$fixed_terms),
    tstart = X_Y$Y[1, ],
    tstop = X_Y$Y[2, ],
    Q_0 = Q_0,

    Q = Q,
    a_0 = a_0,

    Q_tilde = diag(1e-2, n_vars + 1),
    risk_obj = risk_set,
    F = diag(1, n_vars + 1),
    n_max = 10,
    order = 1,
    n_threads = 1,
    N_fw_n_bw = 500,
    N_smooth = 1000,
    N_first = 10000,
    forward_backward_ESS_threshold = NULL,
    debug = 0,
    method = "PF",
    smoother ="Fearnhead_O_N")

  test_func <- function(x){
    cloud <- bquote(.(substitute(x)))
    n <- length(x)
    lapply(1:n, function(i){
      bquote(expect_equal(sum(.(cloud)[[.(i)]]$weights), 1))
    })
  }

  get_test_expr <- function(fit_quote, test_file_name){
    bquote({
      result <- .(fit_quote)

      expect_false(all(old_seed == .Random.seed))

      #####
      # Test versus previous computed values

      # save_to_test(result, file_name = .(test_file_name))
      expect_equal(result, read_to_test(.(test_file_name)), tolerance = 1.49e-08)

      #####
      # Test that weights sum to one
      lapply(test_func(result$forward_clouds), eval)
      lapply(test_func(result$backward_clouds), eval)
      lapply(test_func(result$smoothed_clouds), eval)

      #####
      # Test that multithreaded version gives the same
      old_args <- args
      args$n_threads <- 7

      result_multi <- .(fit_quote)
      expect_equal(result_multi, result)

      args <- old_args
    })
  }

  #####
  # Simple PF
  eval(get_test_expr(
    quote({
      set.seed(30302129)
      old_seed <- .Random.seed
      args$method <- "bootstrap_filter"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/bootstrap_filter"),

    envir = environment())

  #####
  # Normal approximation in importance density with mean from previous cloud
  eval(get_test_expr(
    quote({
      set.seed(30302129)
      old_seed <- .Random.seed
      args$method <- "PF_normal_approx_w_cloud_mean"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/PF_normal_approx_w_cloud_mean"),

    envir = environment())

  #####
  # Normal approximation in AUX filter with mean from previous cloud
  eval(get_test_expr(
    quote({
      set.seed(30302129)
      old_seed <- .Random.seed
      args$method <- "AUX_normal_approx_w_cloud_mean"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/AUX_normal_approx_w_cloud_mean"),

    envir = environment())

  #####
  # Normal approximation in importance density with mean from parent and child particle
  eval(get_test_expr(
    quote({
      set.seed(30302129)
      old_seed <- .Random.seed
      args$method <- "PF_normal_approx_w_particles"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/PF_normal_approx_w_particles"),

    envir = environment())

  #####
  # Normal approximation in AUX filter with mean from parent and child particle
  eval(get_test_expr(
    quote({
      set.seed(30302129)
      old_seed <- .Random.seed
      args$method <- "AUX_normal_approx_w_particles"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/AUX_normal_approx_w_particles"),

    envir = environment())

  #####
  # Test O(n^2) method from Brier et al
  eval(get_test_expr(
    quote({
      set.seed(30302129)
      old_seed <- .Random.seed
      args$method <- "AUX_normal_approx_w_particles"
      args$smoother <- "Brier_O_N_square"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/Brier_O_N_square_AUX_normal_approx_w_particles"),

    envir = environment())

  # # #TODO: clean up
  PF_effective_sample_size <- asNamespace("dynamichazard")$PF_effective_sample_size
  plot(result, type = "smoothed_clouds")
  plot(result, type = "backward_clouds", qlvls = c(), lty = 2, add = TRUE)
  plot(result, type = "forward_clouds", qlvls = c(), lty = 3, add = TRUE)
  abline(h = sims$betas[1, ])
  (tmp <- PF_effective_sample_size(result))
  # (tmp2 <- PF_effective_sample_size(read_to_test("local_tests/AUX_normal_approx_w_particles")))
  # for(i in seq_along(tmp)){
  #   tmp[[i]] <- (tmp[[i]] - tmp2[[i]]) / tmp2[[i]]
  # }
  # tmp
  # ddfit <- ddhazard(
  #   Surv(tstart, tstop, event) ~ . - id,
  #   data = sims$res,
  #   max_T = 10,
  #   by = 1,
  #   id = sims$res$id,
  #   Q_0 = diag(1, 3),
  #   Q = diag(1e-1, 3),
  #   a_0 = sims$betas[1, ],
  #   control = list(NR_eps = 1e-5))
  #
  # sapply(result, function(x){
  #   ws <- lapply(x, "[[", "weights")
  #   sapply(ws, function(z) 1 / sum(z^2))
  # })
  #
  # matplot(0:10, sims$betas, lty = 1, type = "l", ylim = c(-5, 5), xlim = c(0, 11))
  # for(i in 1:3){
  #   state_est <- t(sapply(result[[i]], function(row){
  #     colSums(t(row$states) * drop(row$weights))
  #   }))
  #
  #   idx <- switch(i, "1" = 0:10, "2" = 1:11, "3" = 1:10)
  #   matplot(idx, state_est, lty = i + 1, type = "l", add = TRUE)
  #   matplot(idx, state_est, lty = i + 1, type = "p", add = TRUE, pch = 16 + i)
  # }
  # matplot(0:10, ddfit$state_vecs, lty = 1, col = "blue", type = "l", add = TRUE)
  # # #
  # sapply(lapply(result$backward_clouds, "[[", "parent_idx"),
  #        function(n) sort(xtabs(~ n), decreasing = TRUE)[1:10])
  # sapply(lapply(result$backward_clouds, "[[", "parent_idx"),
  #        function(n) sort(xtabs(~ n), decreasing = TRUE)[1:10])
  #
  # sapply(lapply(result$backward_clouds, "[[", "weights"),
  #        function(x) sort(x, decreasing = TRUE)[1:10])
  # ord <- sapply(lapply(result$backward_clouds, "[[", "weights"),
  #             function(x) order(x, decreasing = TRUE)[1:11])
  # #
  # result$smoothed_clouds[[10]]$states[, ord[, 10]]
  # result$smoothed_clouds[[8]]$states[, 5740]
  #
  # result <- result$smoothed_clouds
})

test_that("Import and export PF cloud from Rcpp gives the same", {
  skip_on_cran()

  PF_cloud_to_rcpp_and_back <- asNamespace("dynamichazard")$PF_cloud_to_rcpp_and_back

  #####
  # Without transition_likelihoods
  cloud_example <- read_to_test("local_tests/cloud_example_no_transition_likelihoods")
  expect_equal( # quick sanity check
    length(cloud_example$transition_likelihoods), 0)

  result <- PF_cloud_to_rcpp_and_back(cloud_example)
  for(i in seq_along(cloud_example)){
    expect_equal(cloud_example[i], result[i])
  }

  #####
  # With transition_likelihoods
  cloud_example <- read_to_test("local_tests/cloud_example_with_transition_likelihoods")
  expect_gt( # quick sanity check
    length(cloud_example$transition_likelihoods), 0)

  result <- PF_cloud_to_rcpp_and_back(cloud_example)
  for(i in seq_along(cloud_example)){
    expect_equal(cloud_example[i], result[i])
  }
})

test_that("compute_summary_stats gives previous results", {
  skip_on_cran()

  compute_summary_stats <- asNamespace("dynamichazard")$compute_summary_stats
  cloud_example <- read_to_test("local_tests/cloud_example_no_transition_likelihoods")

  sum_stats <- compute_summary_stats(cloud_example)

  # save_to_test(sum_stats, "local_tests/compute_summary_stats")
  expect_equal(sum_stats, read_to_test("local_tests/compute_summary_stats"), tolerance = 1.49e-08)

  expect_true(FALSE) # test with O(N^2 method)
})

test_that("PF_EM stops with correct error messages due to wrong or missing arguments", {
  args <- list(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2))

  expect_error(
    do.call(PF_EM, args), paste0(
      "Please supply the number of particle for ", sQuote("control\\$N_first")))

  args$control <- list(N_first = 1e3)
  expect_error(
    do.call(PF_EM, args), paste0(
      "Please supply the number of particle for ", sQuote("control\\$N_fw_n_bw")))

  args$control <- c(args$control, list(N_fw_n_bw = 1e3))
  expect_error(
    do.call(PF_EM, args), paste0(
      "Please supply the number of particle for ", sQuote("control\\$N_smooth")))
})

test_that("PF_EM gives previous results on head neck data set", {
  skip_on_cran()

  set.seed(98612415)
  result <- suppressWarnings( # Supressed as there is a warning about not converging
    PF_EM(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
      control = list(N_fw_n_bw = 200, N_smooth = 1e3, N_first = 2e3,
                     n_max = 5, n_threads = 7),
      max_T = 45))

  #####
  # Test that result are reproducable
  r2 <- result$call
  r2[["control"]] <- c(eval(r2$control, environment()),
                       list(seed = result$seed))
  r2  <- suppressWarnings(eval(r2, environment()))

  result$call <- NULL
  r2$call <- NULL
  expect_equal(result, r2)

  #####
  # Test versus previous computed values
  result <- result[c("a_0", "Q", "clouds", "summary_stats", "log_likes",
                     "n_iter", "effective_sample_size", "seed")]

  # save_to_test(result, "local_tests/PF_head_neck")
  expect_equal(result, read_to_test("local_tests/PF_head_neck"), tolerance = 1.49e-08)
})

test_that("transition_likelihoods are correct and weights sum to one and passing back and forward works",
          expect_true(FALSE))
