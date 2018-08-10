context("Running test_PF")

test_that("dmvnrm_log_test gives correct likelihood", {
  skip_if_not_installed("mvtnorm")

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

# Function to compute cloud means
get_means <- function(result){
  means <- lapply(
    result[c("forward_clouds", "backward_clouds", "smoothed_clouds")],
    function(clouds)
      do.call(rbind, sapply(clouds, function(row){
        colSums(t(row$states) * drop(row$weights))
      }, simplify = FALSE)))
}

test_that("PF_smooth gives same results", {
  skip_on_cran()

  n_vars <- 2
  set.seed(78095324)
  sims <- test_sim_func_logit(
    n_series = 250, n_vars = n_vars, t_0 = 0, t_max = 10,
    x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
    intercept_start = -3, sds = c(.1, rep(.5, n_vars)))

  frm <- Surv(tstart, tstop, event) ~ . - id
  X_Y = get_design_matrix(frm, sims$res)
  risk_set <- get_risk_obj(
    Y = X_Y$Y, by = 1, max_T = 10, id = sims$res$id,
    is_for_discrete_model = TRUE)

  # sum(sims$res$event)
  # matplot(sims$beta, type = "l", lty = 1)
  # xtabs(~ sims$res$tstop[sims$res$event == 1])

  Q <- diag(sqrt(.33), 3)
  Q_0 <- diag(.5, 3)
  a_0 <- sims$betas[1, ]

  args <- list(
    n_fixed_terms_in_state_vec = 0,
    X = t(X_Y$X), fixed_terms = t(X_Y$fixed_terms),
    tstart = X_Y$Y[, 1], tstop = X_Y$Y[, 2], R = diag(1, ncol(Q)), Q_0 = Q_0,
    fixed_parems = numeric(), Q = Q, a_0 = a_0, Q_tilde = diag(1e-2, n_vars + 1),
    risk_obj = risk_set, F = diag(1, n_vars + 1), n_max = 10, n_threads = 1,
    N_fw_n_bw = 20, N_smooth = 100, N_first = 100,
    forward_backward_ESS_threshold = NULL, debug = 0, method = "PF",
    smoother = "Fearnhead_O_N", model = "logit", type = "RW")

  fw_args <- list(
    x = frm, N_fw = args$N_fw_n_bw, N_first = args$N_first, data = sims$res,
    by = 1, fixed_effects = numeric(), control = PF_control(
      N_fw_n_bw = 1, N_smooth = 1, N_first = 1, method = args$method,
      Q_tilde = args$Q_tilde, n_threads = 1),
    max_T = 10, Fmat = diag(1, 3), trace = 0, id = sims$res$id)
  fw_args <- c(fw_args, args[
    intersect(names(args), names(formals(PF_forward_filter.formula)))])

  test_func <- function(test_file_name, update = FALSE){
    get_func <- quote(
      function(x){
        cloud <- bquote(.(substitute(x)))
        n <- length(x)
        lapply(1:n, function(i){
          bquote(expect_equal(sum(.(cloud)[[.(i)]]$weights), 1))
        })
      })

    cl <- list(quote(PF_smooth))
    cl[names(args)] <- lapply(names(args), function(x)
      substitute(args$z, list(z = as.symbol(x))))
    cl <- as.call(cl)

    cl_fw <- list(quote(PF_forward_filter))
    cl_fw[names(fw_args)] <- lapply(names(fw_args), function(x)
      substitute(fw_args$z, list(z = as.symbol(x))))
    cl_fw <- as.call(cl_fw)

    q <- bquote({
      set.seed(30302129)
      old_seed <- .Random.seed

      result <- .(cl)

      expect_false(all(old_seed == .Random.seed))

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
          colSums(result$transition_likelihoods[[2]]$weights)
        })

        for(i in seq_along(ws))
          eval(substitute(
            expect_equal(ws[[i]], rep(1, length(ws[[i]]))), list(i = i)))
      }

      #####
      # We should get the same with the forward filter
      set.seed(30302129)
      fw_res <- .(cl_fw)

      expect_equal(fw_res$forward_clouds, result$forward_clouds)

      #####
      # Changing the seed changes the result
      set.seed(1)
      result_new_seed <- .(cl)

      expect_false(isTRUE(all.equal(result, result_new_seed)))

      #####
      # Test that multithreaded version gives the same
      set.seed(30302129)
      old_args <- args
      args$n_threads <- max(parallel::detectCores(logical = FALSE), 1)

      result_multi <- .(cl)
      expect_equal(result_multi, result)

      args <- old_args

      #####
      # Test versus previous computed values

      # Compute clouds means to test against
      if(.(update)){
        sm <- get_means(result)$smoothed_clouds
        matplot(sm, ylim = range(sm, sims$betas), type = "l", lty = 1)
        matplot(sims$betas, lty = 2, type = "l", add = TRUE)
      }

      .file <- paste0(.(test_file_name), "_cloud_means.RDS")
      eval(substitute(
        expect_known_value(
          get_means(result), .file, tolerance = 1.49e-08, update = .(update)),
        list(.file = .file)), envir = environment())

      # This one may be skipped as the test file is large-ish
      .file <- paste0("local_tests/", .(test_file_name), ".RDS")
      eval(substitute(
        test_if_file_exists(
          .file,
          expect_known_value(
            result, .file, tolerance = 1.49e-08, update = .(update))),
        list(.file = .file)), envir = environment())
    })

    eval(q, envir = parent.frame())
  }

  #####
  # Simple PF
  fw_args$control$method <- args$method <- "bootstrap_filter"
  test_func(test_file_name = "bootstrap_filter")

  #####
  # Normal approximation in importance density with mean from previous cloud
  fw_args$control$method <- args$method <- "PF_normal_approx_w_cloud_mean"
  test_func(test_file_name = "PF_normal_approx_w_cloud_mean")

  #####
  # Normal approximation in AUX filter with mean from previous cloud
  fw_args$control$method <- args$method <- "AUX_normal_approx_w_cloud_mean"
  test_func(test_file_name = "AUX_normal_approx_w_cloud_mean")

  #####
  # Normal approximation in importance density with mean from parent and child particle
  fw_args$control$method <- args$method <- "PF_normal_approx_w_particles"
  test_func(test_file_name = "PF_normal_approx_w_particles")

  #####
  # Normal approximation in AUX filter with mean from parent and child particle
  fw_args$control$method <- args$method <- "AUX_normal_approx_w_particles"
  test_func(test_file_name = "AUX_normal_approx_w_particles")

  #####
  # Test O(n^2) method from Brier et al
  fw_args$control$method <- args$method <- "AUX_normal_approx_w_particles"
  args$smoother <- "Brier_O_N_square"
  test_func(test_file_name = "Brier_AUX_normal_w_particles")
})

test_that("Import and export PF cloud from Rcpp gives the same", {
  skip_on_cran()

  #####
  # Without transition_likelihoods
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

  test_if_file_exists(
    "local_tests/cloud_example_no_transition_likelihoods",
    {
      cloud_example <- read_to_test("local_tests/cloud_example_no_transition_likelihoods")
      test_func(FALSE)
    })

  #####
  # With transition_likelihoods
  test_if_file_exists(
    "local_tests/cloud_example_with_transition_likelihoods",
    {
      cloud_example <- read_to_test("local_tests/cloud_example_with_transition_likelihoods")
      test_func(TRUE)
    })
})

test_that("PF_EM stops with correct error messages due to wrong or missing arguments", {
  args <- list(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer, order = 2, model = "something_not_implemented",
    by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2))

  expect_error(
    do.call(PF_EM, args), paste0(
      sQuote("order"), " not equal to 1 is not supported"))

  args$order <- 1
  expect_error(
    do.call(PF_EM, args), paste0(
      sQuote("model"), " is not supported"))

  args$model <- "logit"
  expect_error(
    do.call(PF_EM, args), paste0(
      "Please supply the number of particle in ",
      sQuote("PF_control(N_first)")), fixed = TRUE)
})

test_that("PF_EM gives previous results on head neck data set", {
  skip_on_cran()

  args <- list(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
    control = list(
      N_fw_n_bw = 200, N_smooth = 1e3, N_first = 2e3, n_max = 5,
      n_threads =  max(parallel::detectCores(logical = FALSE), 1)),
    max_T = 45)

  test_func <- function(
    smoother, file_name, do_update = FALSE, do_print_fname = FALSE){
    f1 <- paste0("local_tests/", file_name)
    f2 <- paste0(file_name, "_cloud_means")

    q <- bquote({
      these_args <- args
      these_args$control <- c(these_args$control, smoother = .(smoother))

      set.seed(98612415)
      # Supressed as there is a warning about not converging
      result <- suppressWarnings(do.call(PF_EM, these_args))

      #####
      # Test that result are reproducable
      r2 <- result$call
      assign(".Random.seed", result$seed, envir = .GlobalEnv)
      r2  <- suppressWarnings(eval(r2, environment()))

      result$call <- NULL
      r2$call <- NULL
      expect_equal(result, r2)

      #####
      # Test versus previous computed values
      result <- result[c("a_0", "Q", "clouds", "EM_ests", "log_likes",
                         "n_iter", "effective_sample_size", "seed")]

      if(.(do_print_fname)){
        func <- function(x)
          cat("tmp <- readRDS('previous_results/", x, ".RDS')\n", sep = "")
        func(.(f1))
        func(.(f2))
      }
      if(.(do_update))
        save_to_test(result, .(f1))

      test_if_file_exists(
        .(f1),
        expect_equal(result, read_to_test(.(f1)), tolerance = 1.49e-08))

      # Compute clouds means to test against
      if(.(do_update))
        save_to_test(get_means(result$clouds), file_name = .(f2))
      expect_equal(get_means(result$clouds), read_to_test(.(f2)),
                   tolerance = 1.49e-08)
    })

    eval(q, envir = parent.frame())
  }

  test_func(smoother = "Fearnhead_O_N", "PF_head_neck")
  test_func(smoother = "Brier_O_N_square", "PF_head_neck_w_Brier_method")
})

test_that("PF_EM gives previous results on head neck data set with fixed effects and the logit model", {
  skip_on_cran()
  skip_if(!dir.exists("previous_results/local_tests"))

  set.seed(61364778)
  pp_fit <- suppressWarnings(
    # we take to few iterations so there will be a warning
    PF_EM(
      formula = survival::Surv(start, stop, event) ~ ddFixed(group),
      data = head_neck_cancer,
      by = 1, Q_0 = 1, Q = 0.01,
      control = PF_control(
        N_fw_n_bw = 200, N_smooth = 2e3, N_first = 2e3,
        n_max = 3,
        method = "AUX_normal_approx_w_cloud_mean",
        n_threads = max(parallel::detectCores(logical = FALSE), 1)),
      max_T = 30))

  # tmp <- readRDS("previous_results/local_tests/pf_logit_w_fixed.RDS")
  expect_known_value(pp_fit[!names(pp_fit) %in% c("clouds", "call")],
                     "local_tests/pf_logit_w_fixed.RDS")
})

test_that("compute_summary_stats_first_o_RW gives previous results", {
  skip_on_cran()

  Q_0 <- diag(1, 2)
  Q <- diag(0.1, 2)
  a_0 <- c(-3.6, 0)
  R <- diag(1, 2)

  test_func <- function(data_file, test_file, do_update = FALSE){
    q <- bquote({
      test_if_file_exists(
        .(data_file),
        {
          cloud_example <- read_to_test(.(data_file))$clouds

          #####
          # Test that multithreaded version gives the same
          sum_stats <- compute_summary_stats_first_o_RW(
            cloud_example, 1, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE)
          s2 <- compute_summary_stats_first_o_RW(
            cloud_example, 4, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE)

          expect_equal(sum_stats, s2)

          #####
          # Test versus previous results
          if(.(do_update))
            save_to_test(sum_stats, .(test_file))
          expect_equal(sum_stats, read_to_test(.(test_file)), tolerance = 1.49e-08)
      })
    })

    eval(q, envir = parent.frame())
  }

  test_func("local_tests/PF_head_neck", "local_tests/compute_summary_stats")

  ######
  # Test with method from Brier et al (2010)
  test_func("local_tests/PF_head_neck_w_Brier_method",
            "local_tests/compute_summary_stats_w_Brier_method")
})

test_that("´est_params_dens´ gives the same as a R version", {
  skip_on_cran()
  fit <- suppressWarnings(PF_EM(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
    control = list(
      N_fw_n_bw = 200, N_smooth = 1e3, N_first = 2e3, n_max = 1),
    max_T = 45, type = "VAR"))

  . <- function(bw, this){
    ws <- sqrt(this$weights)
    good <- which(ws > 1e-8)
    ws <- ws[good]

    X_new <- t(bw$states[, this$parent_idx[good], drop = FALSE])
    Y_new <- t(this$states[, good, drop = FALSE])

    list(X = X_new, Y = Y_new, sqrt_ws = ws)
  }

  X <- Y <- matrix(ncol = 2)[-1, , drop = FALSE]
  sqrt_ws <- numeric()
  fw_clouds <- fit$clouds$forward_clouds
  clouds <- fit$clouds$smoothed_clouds

  for(i in 2:length(clouds)){
    out <- .(fw_clouds[[i]], clouds[[i]])

    sqrt_ws <- c(sqrt_ws, out$sqrt_ws)
    X <- rbind(X, out$X)
    Y <- rbind(Y, out$Y)
  }

  wX <- X * sqrt_ws
  wY <- Y * sqrt_ws
  lm_fit <- lm(Y ~ X - 1, weights = sqrt_ws^2)
  expect_equal(lm_fit$coefficients, t(fit$F), check.attributes = FALSE)

  qr_obj <- qr(wX)
  qty <- qr.qty(qr_obj, wY)[1:ncol(wX), ]
  Q <- (crossprod(wY) - crossprod(qty)) / (length(clouds) - 1)
  expect_equal(Q, fit$Q, check.attributes = FALSE)

  Q <- crossprod(wY - wX %*% lm_fit$coefficients) / (length(clouds) - 1)
  expect_equal(Q, fit$Q, check.attributes = FALSE)
})

test_that("fixed effect estimation gives the same as an R implementation", {
  skip_on_cran()

  # specific function for this example. Only works when 1 fixed and 1 varying
  # parameter
  R_est_func <- function(pf_res, fam, start){
    cl <- pf_res$call
    cl <- cl[c(1, match(
      c("formula", "data", "max_T", "id", "trace", "order", "by", "model"),
      names(cl), 0))]
    cl$trace <- 0
    cl[[1]] <- quote(dynamichazard:::.get_PF_static_args)
    static_args <- eval(cl)

    cls <- pf_res$clouds$smoothed_clouds
    d <- length(cls)
    tstop  <- static_args$tstop
    tstart <- static_args$tstart

    y <- offs <- Z <- ws <- NULL
    for(i in 1:d){
      r_set <- static_args$risk_obj$risk_sets[[i]]
      X_i <- static_args$X[, r_set]
      Z_i <- static_args$fixed_terms[, r_set]
      Y_i <- static_args$risk_obj$is_event_in[r_set] == (i - 1L)

      cl <- cls[[i]]
      good <- which(drop(cl$weights > 1e-7))
      ws_i <- cl$weights[good]
      sta <- cl$states[good]

      off <- c(sapply(sta, function(s) s * X_i))
      if(fam == "exponential")
        off <- off + log(pmin(tstop[r_set], i) - pmax(tstart[r_set], i - 1))

      ws_i <- c(sapply(ws_i, rep, times = length(r_set)))
      Z_i <- rep(Z_i, length(good))
      Y_i <- rep(Y_i, length(good))

      ws <- c(ws, ws_i)
      Z <- c(Z, Z_i)
      offs <- c(offs, off)
      y <- c(y, Y_i)
    }

    fam <- if(fam == "exponential") poisson else binomial
    suppressWarnings(
      # avoid `non-integer #successes in a binomial glm!`
      # and the non-covergence warning
      glm(y ~ Z - 1, family = fam, offset = offs, weights = ws,
          control = glm.control(maxit = 1), start = start))
  }

  fam <- "logit"
  set.seed(19724898)
  pf_fit <- suppressWarnings(PF_EM(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer, id = head_neck_cancer$id, model = fam,
    by = 1, Q_0 = 1, Q = .1,
    control = list(
      N_fw_n_bw = 100, N_smooth = 200, N_first = 1000, n_max = 1),
    max_T = 45, type = "RW"))

  # start is found by running the above with `trace = 1`
  fit <- R_est_func(pf_fit, fam, start = 0.6576882362117193)
  expect_equal(
    fit$coefficients, pf_fit$fixed_effects, check.attributes = FALSE)

  fam <- "exponential"
  pf_fit <- suppressWarnings(PF_EM(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer, id = head_neck_cancer$id, model = fam,
    by = 1, Q_0 = 1, Q = .1,
    control = list(
      N_fw_n_bw = 100, N_smooth = 200, N_first = 1000, n_max = 1),
    max_T = 45, type = "RW"))
  fit <- R_est_func(pf_fit, fam, start = 0.6255092546830791)
  expect_equal(
    fit$coefficients, pf_fit$fixed_effects, check.attributes = FALSE)
})
