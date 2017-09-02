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
    N_fw_n_bw = 20,
    N_smooth = 100,
    N_first = 100,
    forward_backward_ESS_threshold = NULL,
    debug = 0,
    method = "PF",
    smoother = "Fearnhead_O_N")

  test_func <- function(fit_quote, test_file_name){
    get_func <- quote(
      function(x){
        cloud <- bquote(.(substitute(x)))
        n <- length(x)
        lapply(1:n, function(i){
          bquote(expect_equal(sum(.(cloud)[[.(i)]]$weights), 1))
        })
      })

    q <- bquote({
      set.seed(30302129)
      old_seed <- .Random.seed
      result <- .(fit_quote)
      expect_false(all(old_seed == .Random.seed))

      #####
      # Test versus previous computed values

      # save_to_test(result, file_name = .(test_file_name))
      expect_equal(result, read_to_test(.(test_file_name)), tolerance = 1.49e-08)

      #####
      # Test that weights sum to one
      func <- .(get_func)

      lapply(func(result$forward_clouds), eval, envir = environment())
      lapply(func(result$backward_clouds), eval, envir = environment())
      lapply(func(result$smoothed_clouds), eval, envir = environment())

      #####
      # Check that weights in transition_likelihoods sum to one
      if(length(result$transition_likelihoods) > 0){
        ws <- lapply(result$transition_likelihoods, function(x){
          sapply(lapply(x, "[[", "weights"), sum)
        })

        for(i in seq_along(ws))
          eval(substitute(
            expect_equal(ws[[i]], rep(1, length(ws[[i]]))), list(i = i)))
      }

      #####
      # Test that multithreaded version gives the same
      old_args <- args
      args$n_threads <- 7

      result_multi <- .(fit_quote)
      expect_equal(result_multi, result)

      args <- old_args
    })

    eval(q, envir = parent.frame())
  }

  #####
  # Simple PF
  test_func(
    quote({
      args$method <- "bootstrap_filter"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/bootstrap_filter")

  #####
  # Normal approximation in importance density with mean from previous cloud
  test_func(
    quote({
      args$method <- "PF_normal_approx_w_cloud_mean"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/PF_normal_approx_w_cloud_mean")

  #####
  # Normal approximation in AUX filter with mean from previous cloud
  test_func(
    quote({
      args$method <- "AUX_normal_approx_w_cloud_mean"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/AUX_normal_approx_w_cloud_mean")

  #####
  # Normal approximation in importance density with mean from parent and child particle
  test_func(
    quote({
      args$method <- "PF_normal_approx_w_particles"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/PF_normal_approx_w_particles")

  #####
  # Normal approximation in AUX filter with mean from parent and child particle
  test_func(
    quote({
      args$method <- "AUX_normal_approx_w_particles"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/AUX_normal_approx_w_particles")

  #####
  # Test O(n^2) method from Brier et al
  test_func(
    quote({
      args$method <- "AUX_normal_approx_w_particles"
      args$smoother <- "Brier_O_N_square"
      result <- do.call(PF_smooth, args)
      result
    }),
    test_file_name = "local_tests/Brier_O_N_square_AUX_normal_approx_w_particles")
})

test_that("Import and export PF cloud from Rcpp gives the same", {
  skip_on_cran()

  PF_cloud_to_rcpp_and_back <- asNamespace("dynamichazard")$PF_cloud_to_rcpp_and_back

  #####
  # Without transition_likelihoods
  cloud_example <- read_to_test("local_tests/cloud_example_no_transition_likelihoods")

  test_func <- function(has_transition){
    q <- bquote({
      .(if(has_transition)
        as.symbol("expect_gt")
        else as.symbol("expect_equal"))( # quick sanity check
        length(cloud_example$transition_likelihoods), 0)

      result <- PF_cloud_to_rcpp_and_back(cloud_example)
      for(i in seq_along(cloud_example)){
        eval(substitute(
          expect_equal(cloud_example[[i]], result[[i]]),
          list(i = i)), envir = environment())
      }
    })

    eval(q, envir = parent.frame())
  }

  test_func(FALSE)

  #####
  # With transition_likelihoods
  cloud_example <- read_to_test("local_tests/cloud_example_with_transition_likelihoods")

  test_func(TRUE)
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

  args <- list(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
    control = list(N_fw_n_bw = 200, N_smooth = 1e3, N_first = 2e3,
                   n_max = 5, n_threads = 7),
    max_T = 45)

  test_func <- function(smoother, file_name){
    q <- bquote({
      these_args <- args
      these_args$control <- c(these_args$control, smoother = .(smoother))

      set.seed(98612415)
      result <- suppressWarnings( # Supressed as there is a warning about not converging
        do.call(PF_EM, these_args))

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

      # save_to_test(result, .(file_name))
      expect_equal(result, read_to_test(.(file_name)), tolerance = 1.49e-08)
    })

    eval(q, envir = parent.frame())
  }

  test_func(smoother = "Fearnhead_O_N", "local_tests/PF_head_neck")
  test_func(smoother = "Brier_O_N_square", "local_tests/PF_head_neck_w_Brier_method")
})

test_that("compute_summary_stats gives previous results", {
  skip_on_cran()

  compute_summary_stats <- asNamespace("dynamichazard")$compute_summary_stats

  Q_0 <- diag(1, 2)
  Q <- diag(0.1, 2)
  a_0 <- c(-3.6, 0)

  test_func <- function(data_file, test_file){
    q <- bquote({
      cloud_example <- read_to_test(.(data_file))$clouds

      #####
      # Test that multithreaded version gives the same
      sum_stats <- compute_summary_stats(cloud_example, 1, a_0 = a_0, Q = Q, Q_0 = Q_0)
      s2 <- compute_summary_stats(cloud_example, 4, a_0 = a_0, Q = Q, Q_0 = Q_0)

      expect_equal(sum_stats, s2)

      #####
      # Test versus previous results
      # save_to_test(sum_stats, .(test_file))
      expect_equal(sum_stats, read_to_test(.(test_file)), tolerance = 1.49e-08)
    })

    eval(q, envir = parent.frame())
  }

  test_func("local_tests/PF_head_neck", "local_tests/compute_summary_stats")

  ######
  # Test with method from Brier et al (2010)
  test_func("local_tests/PF_head_neck_w_Brier_method",
            "local_tests/compute_summary_stats_w_Brier_method")
})

test_that("help page example runs and gives previous computed results", {
  skip_on_cran()

  .lung <- lung[!is.na(lung$ph.ecog), ]
  pf_fit <- PF_EM(
    Surv(time, status == 2) ~ ph.ecog + age,
    data = .lung, by = 50, id = 1:nrow(.lung),
    Q_0 = diag(1, 3), Q = diag(1, 3),
    max_T = 800,
    control = list(
      N_fw_n_bw = 500,
      N_first = 2500,
      N_smooth = 2500,
      n_max = 25,
      n_threads = parallel::detectCores()),
    trace = 1)

  # save_to_test(pf_fit[names(pf_fit) != "call"], "local_tests/survival_lung_example")
  expect_equal(pf_fit[names(pf_fit) != "call"], read_to_test("local_tests/survival_lung_example"), tolerance = 1.49e-08)

  expect_no_error(plot(pf_fit, cov_index = 1))
  expect_no_error(plot(pf_fit, cov_index = 2))
  expect_no_error(plot(pf_fit, cov_index = 3))
  expect_no_error(plot(pf_fit$log_likes))
})
