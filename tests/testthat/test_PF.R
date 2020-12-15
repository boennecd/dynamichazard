context("Running test_PF")

test_that("dmvnrm_log_test gives correct likelihood", {
  skip_if_not_installed("mvtnorm")
  library(mvtnorm)

  set.seed(44518266)
  for(i in 1:10){
    n <- 5
    mu <- rnorm(n)
    A <- matrix(runif(n^2)*2-1, ncol=n)
    sigma <- t(A) %*% A
    sigma_chol_inv <- solve(chol(sigma))

    for(i in 1:10){
      x <- rnorm(n)
      expect_equal(
        dmvnorm(x, mu, sigma, log = TRUE),
        dmvnrm_log_test(x, mu, sigma_chol_inv))
    }
  }
})

test_that("dmvtrm_log_test gives correct likelihood", {
  skip_if_not_installed("mvtnorm")
  library(mvtnorm)

  # dfunc <- function(x, mu, sigma, nu){
  #   p <- length(x)
  #   constant <-
  #     lgamma((nu + p) / 2) - lgamma(nu / 2) -
  #     log(nu * pi) * (p / 2) - log(det(sigma)) /2
  #   delta <- x - mu
  #   cp <- crossprod(delta, solve(sigma, delta))
  #
  #   drop(constant - (nu + p) / 2 * log1p(cp / nu))
  # }

  set.seed(44518266)
  for(i in 1:10){
    n <- 5
    mu <- rnorm(n)
    A <- matrix(runif(n^2)*2-1, ncol=n)
    sigma <- t(A) %*% A
    sigma_chol_inv <- solve(chol(sigma))
    nu <- i + 2L

    for(i in 1:10){
      x <- rnorm(n)
      expect_equal(
        dmvt(x = x, delta = mu, sigma = sigma, df = nu, log = TRUE),
        dmvtrm_log_test(x, mu, sigma_chol_inv = sigma_chol_inv, nu = nu))
      # expect_equal(
      #   dfunc(x = x, mu = mu, sigma = sigma, nu = nu),
      #   dmvtrm_log_test(x, mu, sigma_chol_inv = sigma_chol_inv, nu = nu))
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
    fixed_params = numeric(), Q = Q, a_0 = a_0, Q_tilde = diag(1e-2, n_vars + 1),
    risk_obj = risk_set, F = diag(1, n_vars + 1), n_max = 10, n_threads = 1,
    N_fw_n_bw = 20, N_smooth = 100, N_first = 100, N_smooth_final = 100,
    forward_backward_ESS_threshold = NULL, debug = 0, ftol_rel = 1e-8,
    method = "bootstrap_filter", covar_fac = -1,
    smoother = "Fearnhead_O_N", model = "logit", type = "RW", nu = 0L)

  fw_args <- list(
    x = frm, N_fw = args$N_fw_n_bw, N_first = args$N_first, data = sims$res,
    by = 1, fixed_effects = numeric(), control = PF_control(
      N_fw_n_bw = 1, N_smooth = 1, N_first = 1, method = args$method,
      Q_tilde = args$Q_tilde, n_threads = 1),
    max_T = 10, trace = 0, id = sims$res$id)
  fw_args <- c(fw_args, args[
    intersect(names(args), names(formals(PF_forward_filter.formula)))])

  test_func <- function(test_file_name, update = do_update_tests){
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

  # .lung <- lung[!is.na(lung$ph.ecog), ]
  # .lung$age <- scale(.lung$age)
  # # fit
  # set.seed(43588155)
  # pf_fit <- PF_EM(
  #   Surv(time, status == 2) ~ ddFixed(ph.ecog) + age,
  #   data = .lung, by = 50, id = 1:nrow(.lung),
  #   Q_0 = diag(1, 2), Q = diag(.5^2, 2),
  #   max_T = 600,
  #   control = PF_control(
  #     N_fw_n_bw = 100, N_first = 200, N_smooth = 100,
  #     n_max = 1, nu = 5L, est_a_0 = FALSE,
  #     smoother = "Fearnhead_O_N", method = "bootstrap_filter",
  #     n_threads = max(parallel::detectCores(logical = FALSE), 1)))
  # saveRDS(
  #   pf_fit$clouds,
  #   "previous_results/local_tests/cloud_example_no_transition_likelihoods.RDS")
  # pf_fit <- PF_EM(
  #   Surv(time, status == 2) ~ ddFixed(ph.ecog) + age,
  #   data = .lung, by = 50, id = 1:nrow(.lung),
  #   Q_0 = diag(1, 2), Q = diag(.5^2, 2),
  #   max_T = 600,
  #   control = PF_control(
  #     N_fw_n_bw = 100, N_first = 200, N_smooth = 100,
  #     n_max = 1, nu = 5L, est_a_0 = FALSE,
  #     smoother = "Brier_O_N_square", method = "bootstrap_filter",
  #     n_threads = max(parallel::detectCores(logical = FALSE), 1)))
  # saveRDS(
  #   pf_fit$clouds,
  #   "previous_results/local_tests/cloud_example_with_transition_likelihoods.RDS")

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
    smoother, file_name, do_update = do_update_tests, do_print_fname = FALSE){
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
  expect_known_value(pp_fit[!names(pp_fit) %in% c("clouds", "call", "control")],
                     "local_tests/pf_logit_w_fixed.RDS")
  expect_equal(
    pp_fit$control$n_threads, max(parallel::detectCores(logical = FALSE), 1))
})

test_that("compute_PF_summary_stats gives previous results", {
  skip_on_cran()

  F. <- diag(1, 2)
  F1 <- matrix(c(.5, 0, .1, .3), 2)
  Q_0 <- diag(1, 2)
  Q <- diag(0.1, 2)
  a_0 <- c(-3.6, 0)
  R <- diag(1, 2)

  test_func <- function(data_file, test_file, do_update = do_update_tests){
    q <- bquote({
      test_if_file_exists(
        .(data_file),
        {
          cloud_example <- read_to_test(.(data_file))$clouds

          #####
          # Test that multithreaded version gives the same
          sum_stats <- compute_PF_summary_stats(
            cloud_example, 1, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE, F = F.)
          s2 <- compute_PF_summary_stats(
            cloud_example, 4, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE, F = F.)

          expect_equal(sum_stats, s2)

          # we should get the same when we use F
          s3 <- compute_PF_summary_stats(
            cloud_example, 4, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE, F = F., do_use_F = TRUE)
          expect_equal(sum_stats, s3)

          # should only compute the latter sets of elements only
          s4 <- compute_PF_summary_stats(
            cloud_example, 4, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE, F = F., do_compute_E_x = FALSE)
          expect_equal(lapply(sum_stats, "[[", "E_x_less_x_less_one_outers"),
                       lapply(s4       , "[[", "E_x_less_x_less_one_outers"))
          expect_true(all(do.call(rbind, lapply(s4, "[[", "E_xs")) == 0))

          # should get the same when we compute with F on our site
          cp <- cloud_example
          for(i in 1:(length(cp$forward_clouds) - 1))
            cp$forward_clouds[[i]]$states <-
              F1 %*% cp$forward_clouds[[i]]$states
          s5 <- compute_PF_summary_stats(
            cloud_example, 4, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE, F = F1, do_use_F = TRUE)
          s6 <- compute_PF_summary_stats(
            cp           , 4, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
            debug = FALSE, F = F., do_use_F = FALSE)
          expect_equal(s5, s6)

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

test_that("'get_cloud_means' and 'get_cloud_quantiles' gives previous results", {
  skip_on_cran()
  skip_if(!dir.exists("previous_results/local_tests"))

  pf_fit <- read_to_test("local_tests/PF_head_neck")
  class(pf_fit) <- "PF_EM"

  #####
  # means
  m1 <- get_cloud_means(pf_fit)
  m2 <- get_cloud_means(pf_fit$clouds)
  expect_equal(m1, m2)
  expect_known_value(m1, file = "PF_mean_2D.RDS")

  m3 <- get_cloud_means(pf_fit, cov_index = 1L)
  expect_equal(dim(m3), c(nrow(m1), 1L))
  expect_equal(m3, m1[, 1L, drop = FALSE])

  m4 <- get_cloud_means(pf_fit, type = "forward_clouds")
  expect_true(!isTRUE(all.equal(m1, m4)))

  m5 <- get_cloud_means(pf_fit, type = "backward_clouds")
  expect_true(!isTRUE(all.equal(m1, m5)))

  #####
  # quantiles
  q1 <- get_cloud_quantiles(pf_fit)
  q2 <- get_cloud_quantiles(pf_fit$clouds)
  expect_equal(q1, q2)
  expect_known_value(q1, file = "PF_quantile_2D.RDS")

  q3 <- get_cloud_quantiles(pf_fit, cov_index = 1L, qlvls = .5)
  expect_equal(dim(q3), c(1L, 1L, dim(q1)[3]))
  expect_equal(q3, q1[2L, 1L, , drop = FALSE])

  q4 <- get_cloud_quantiles(pf_fit, qlvls = .5)
  expect_equal(dim(q4), c(1L, dim(q1)[2:3]))
  expect_equal(q4, q1[2L, , , drop = FALSE])

  q5 <- get_cloud_quantiles(pf_fit, cov_index = 1L)
  expect_equal(dim(q5), c(dim(q1)[1L], 1L, dim(q1)[3L]))
  expect_equal(q5, q1[, 1L, , drop = FALSE])

  q6 <- get_cloud_quantiles(pf_fit, type = "forward_clouds")
  expect_true(!isTRUE(all.equal(q1, q6)))

  q7 <- get_cloud_quantiles(pf_fit, type = "backward_clouds")
  expect_true(!isTRUE(all.equal(q1, q7)))
})

test_that("'get_ancestors' yields the correct result", {
  skip_on_cran()
  skip_if(!dir.exists("previous_results/local_tests"))

  pf_fit <- read_to_test("local_tests/PF_head_neck")
  class(pf_fit) <- "PF_EM"

  cpp_out <- test_get_ancestors(pf_fit$clouds)

  fw <- pf_fit$clouds$forward_clouds
  n_ele <- length(fw)
  correct <- list()
  correct[[n_ele]] <- seq_len(length(fw[[n_ele]]$weights))
  for(i in rev(seq_len(n_ele - 1L)))
    correct[[i]] <- unique(fw[[i + 1L]]$parent_idx[correct[[i + 1L]]])

  for(i in seq_along(correct))
    correct[[i]] <- correct[[i]] - 1L

  expect_equal(correct, cpp_out)
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

  for(i in 1:length(clouds)){
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
  Q <- (crossprod(wY) - crossprod(qty)) / length(clouds)
  expect_equal(Q, fit$Q, check.attributes = FALSE)

  Q <- crossprod(wY - wX %*% lm_fit$coefficients) / length(clouds)
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

test_that("A few iterations with `type = \"VAR\"' yields the same as before", {
  skip_on_cran()

  n_obs     <- 2000L
  n_periods <- 300L

  # Fmat <- matrix(c(.8, 0, 0, .8), 2)
  # Rmat <- diag(1    , 2)
  # Qmat <- diag(.33^2, 2)
  # Q_0  <- get_Q_0(Qmat, Fmat)
  # beta <- c(-6.5, -2)
  #
  # set.seed(54432125)
  # betas <- matrix(nrow = n_periods + 1, ncol = 2)
  # betas[1, ] <- rnorm(2) %*% chol(Q_0)
  # for(i in 1:n_periods + 1)
  #   betas[i, ] <- Fmat %*% betas[i - 1, ] + drop(rnorm(2) %*% chol(Qmat))
  #
  # betas <- t(t(betas) + beta)
  #
  # df <- replicate(n_obs, {
  #   # left-censoring
  #   tstart <- max(0L, sample.int((n_periods - 1L) * 2L, 1) - n_periods + 1L)
  #
  #   # covariates
  #   x <- runif(1, -1, 1)
  #   covars <- c(1, x)
  #
  #   # outcome (stop time and event indicator)
  #   y <- FALSE
  #   for(tstop in (tstart + 1L):n_periods){
  #     fail_time <- rexp(1) / exp(covars %*% betas[tstop + 1L, ])
  #     if(fail_time <= 1){
  #       y <- TRUE
  #       tstop <- tstop - 1L + fail_time
  #       break
  #     }
  #   }
  #
  #   c(tstart = tstart, tstop = tstop, x = x, y = y)
  # })
  # df <- data.frame(t(df))
  # saveRDS(df, "VAR_data.RDS")
  df <- readRDS("VAR_data.RDS")

  # this tend toward the true values. We only take a few iterations though...
  set.seed(30520116)
  pf_Fear <- suppressWarnings(PF_EM(
    Surv(tstart, tstop, y) ~ x + ddFixed(x) + ddFixed_intercept(TRUE), df,
    Q_0 = diag(1, 2), Q = diag(1, 2), Fmat = matrix(c(.1, 0, 0, .1), 2),
    by = 1, type = "VAR", model = "exponential", max_T = n_periods,
    control = PF_control(
      N_fw_n_bw = 50, N_smooth = 100, N_first = 500, eps = .001,
      method = "AUX_normal_approx_w_cloud_mean",
      n_max = 3, smoother = "Fearnhead_O_N",
      Q_tilde = diag(.3^2, 2), n_threads = 4)))

  # old <- readRDS(file.path("previous_results", "PF_VARS.RDS"))
  expect_known_value(pf_Fear[!names(pf_Fear) %in% c("clouds", "terms")],
                     file = "PF_VARS.RDS",
                     # low tol because of a large difference on Solaris.
                     # Likely due to re-sampling being triggered at some point
                     tolerance = 1e-3)
})

test_that("Commutation matrix is correct", {
  A <- structure(c(0.96, 0.27, -0.25, -1.67, 0.63, 0.12, 1.29, 0.31),
                 .Dim = c(4L, 2L))
  expect_equal(as.vector(t(A)), drop(.get_cum_mat(4, 2) %*% as.vector(A)))
  expect_equal(as.vector(A), drop(.get_cum_mat(2, 4) %*% as.vector(t(A))))
})

test_that("PF_EM gives the same with restricted and unrestricted model when we estimate all the parameters", {
  skip_on_cran()

  n_obs     <- 200L
  n_periods <- 100L

  # Fmat <- matrix(c(.9, 0, 0, .9), 2)
  # Rmat <- diag(1    , 2)
  # Qmat <- diag(.33^2, 2)
  # Q_0  <- get_Q_0(Qmat, Fmat)
  # beta <- c(-6.5, -2)
  #
  # set.seed(54432125)
  # betas <- matrix(nrow = n_periods + 1, ncol = 2)
  # betas[1, ] <- rnorm(2) %*% chol(Q_0)
  # for(i in 1:n_periods + 1)
  #   betas[i, ] <- Fmat %*% betas[i - 1, ] + drop(rnorm(2) %*% chol(Qmat))
  #
  # betas <- t(t(betas) + beta)
  #
  # df <- replicate(n_obs, {
  #   # left-censoring
  #   tstart <- max(0L, sample.int((n_periods - 1L) * 2L, 1) - n_periods + 1L)
  #
  #   # covariates
  #   x <- runif(1, -1, 1)
  #   covars <- c(1, x)
  #
  #   # outcome (stop time and event indicator)
  #   y <- FALSE
  #   for(tstop in (tstart + 1L):n_periods){
  #     fail_time <- rexp(1) / exp(covars %*% betas[tstop + 1L, ])
  #     if(fail_time <= 1){
  #       y <- TRUE
  #       tstop <- tstop - 1L + fail_time
  #       break
  #     }
  #   }
  #
  #   c(tstart = tstart, tstop = tstop, x = x, y = y)
  # })
  # df <- data.frame(t(df))
  # saveRDS(df, "VAR_restrict_vs_non_data.RDS")
  df <- readRDS("VAR_restrict_vs_non_data.RDS")

  set.seed(seed <- 30520116)
  pf_non_restrict_Fear <- suppressWarnings(PF_EM(
    Surv(tstart, tstop, y) ~ x + ddFixed(x) + ddFixed_intercept(TRUE), df,
    Q_0 = diag(1, 2), Q = diag(1, 2), Fmat = matrix(c(.1, 0, 0, .1), 2),
    by = 1, type = "VAR", model = "exponential", max_T = n_periods,
    control = PF_control(
      N_fw_n_bw = 50, N_smooth = 100, N_first = 500, eps = .001,
      method = "AUX_normal_approx_w_cloud_mean",
      n_max = 3, smoother = "Fearnhead_O_N",
      Q_tilde = diag(.3^2, 2), n_threads = 4)))

  cl <- pf_non_restrict_Fear$call
  cl[c("Q", "Fmat")] <- NULL
  cl[c("G", "J", "K", "theta", "psi", "phi")] <- list(
    G = diag(2^2), J = diag(2), K = diag(2 * (2 - 1) / 2),
    theta = c(.1, 0, 0, .1), psi = c(0, 0), phi = 0)

  set.seed(seed)
  pf_restrict_Fear <- suppressWarnings(eval(cl))

  expect_equal(
    pf_non_restrict_Fear$F, pf_restrict_Fear$F, tolerance = 1e-6)
  expect_equal(
    pf_non_restrict_Fear$Q, pf_restrict_Fear$Q, tolerance = 1e-5)
})

test_that("type = 'VAR' works with non-zero mean with a single term and gives previous results", {
  skip_on_os("solaris")

  # had somes issues when with an example like this
  set.seed(30520116)
  pf_Fear <- suppressWarnings(PF_EM(
    Surv(start, stop, event) ~ ddFixed(group) + ddFixed_intercept(TRUE),
    head_neck_cancer, Q_0 = 1, Q = .2, Fmat = diag(1e-2, 1),
    by = 1, type = "VAR", model = "logit", max_T = 30,
    control = PF_control(
      N_fw_n_bw = 50, N_smooth = 100, N_first = 500, eps = .001,
      method = "AUX_normal_approx_w_cloud_mean",
      n_max = 2, smoother = "Fearnhead_O_N",
      Q_tilde = diag(.3^2, 1), n_threads = 1)))

  # old <- readRDS(file.path(
  #    "previous_results", "PF_VARS_non_zero_mean_inter.RDS"))
  expect_known_value(pf_Fear[!names(pf_Fear) %in% c("clouds", "terms")],
                     file = "PF_VARS_non_zero_mean_inter.RDS",

                     )

  set.seed(30520116)
  pf_Fear <- suppressWarnings(PF_EM(
    Surv(start, stop, event) ~ group + ddFixed(group) + ddFixed_intercept(),
    head_neck_cancer, Q_0 = 1, Q = .2, Fmat = diag(1e-2, 1),
    by = 1, type = "VAR", model = "logit", max_T = 30,
    control = PF_control(
      N_fw_n_bw = 50, N_smooth = 100, N_first = 500, eps = .001,
      method = "AUX_normal_approx_w_cloud_mean",
      n_max = 2, smoother = "Fearnhead_O_N",
      Q_tilde = diag(.3^2, 1), n_threads = 1)))

  # old <- readRDS(file.path(
  #    "previous_results", "PF_VARS_non_zero_mean_slope.RDS"))
  expect_known_value(pf_Fear[!names(pf_Fear) %in% c("clouds", "terms")],
                     file = "PF_VARS_non_zero_mean_slope.RDS",
                     # low tol because of a large difference on Solaris.
                     # Likely due to re-sampling being triggered at some point
                     tolerance = 1e-3)
})

test_that("Using `n_smooth_final` works as expected and yields previous results", {
  skip_on_cran()

  set.seed(seed <- 56219373)
  .lung <- lung[!is.na(lung$ph.ecog), ]
  .lung$age <- scale(.lung$age)
  f_fit_1 <- suppressWarnings(PF_EM(
    Surv(time, status == 2) ~ ddFixed(ph.ecog) + age,
    data = .lung, by = 50, id = 1:nrow(.lung),
    Q_0 = diag(1, 2), Q = diag(.5^2, 2),
    max_T = 800,
    control = PF_control(
      N_fw_n_bw = 500, N_first = 2500, N_smooth = 5000, N_smooth_final = 1000,
      n_max = 1, eps = .001, Q_tilde = diag(.2^2, 2), est_a_0 = FALSE,
      n_threads = 1)))

  expect_known_value(f_fit_1[!names(f_fit_1) %in% c("clouds", "call", "terms")],
                     file = "n_smooth_final_RW.RDS")

  # # compare with the following after you change `n_max`
  # set.seed(seed)
  # pf_fit <- PF_EM(
  #   Surv(time, status == 2) ~ ddFixed(ph.ecog) + age,
  #   data = .lung, by = 50, id = 1:nrow(.lung),
  #   Q_0 = diag(1, 2), Q = diag(.5^2, 2),
  #   max_T = 800,
  #   control = PF_control(
  #     N_fw_n_bw = 500, N_first = 2500, N_smooth = 5000,
  #     n_max = 50, eps = .001, Q_tilde = diag(.2^2, 2), est_a_0 = FALSE,
  #     n_threads = max(parallel::detectCores(logical = FALSE), 1)), trace = 1)

  # plot(pf_fit$log_likes)
  # points(seq_along(f_fit_1$log_likes), f_fit_1$log_likes, pch = 16)

  set.seed(seed)
  f_fit_2 <- suppressWarnings(PF_EM(
    fixed = Surv(time, status == 2) ~ ph.ecog + age, random = ~ age,
    data = .lung, by = 50, id = 1:nrow(.lung),
    Q_0 = diag(1, 2), Q = diag(.5^2, 2), Fmat = diag(.9, 2),
    max_T = 800, type = "VAR",
    control = PF_control(
      N_fw_n_bw = 500, N_first = 2500, N_smooth = 1000, N_smooth_final = 500,
      n_max = 1, eps = .001, Q_tilde = diag(.1^2, 2),
      n_threads = 1)))

  # # compare with the following after you change `n_max`
  # set.seed(seed)
  # pf_fit_2 <- suppressWarnings(PF_EM(
  #   fixed = Surv(time, status == 2) ~ ph.ecog + age, random = ~ age,
  #   data = .lung, by = 50, id = 1:nrow(.lung),
  #   Q_0 = diag(1, 2), Q = diag(.5^2, 2), Fmat = diag(.9, 2),
  #   max_T = 800, type = "VAR",
  #   control = PF_control(
  #     N_fw_n_bw = 1000, N_first = 2500, N_smooth = 5000,
  #     n_max = 50, eps = .001, Q_tilde = diag(.1^2, 2),
  #     n_threads = max(parallel::detectCores(logical = FALSE), 1)),
  #   trace = 1))
  # plot(pf_fit_2$log_likes)
  # points(seq_along(f_fit_2$log_likes), f_fit_2$log_likes, pch = 16)

  # old <- readRDS(file.path(
  #    "previous_results", "n_smooth_final_VAR.RDS"))
  expect_known_value(
    f_fit_2[!names(f_fit_2) %in% c("clouds", "call", "terms", "fixed",
                                   "random")],
    file = "n_smooth_final_VAR.RDS")
})

test_that("sampling with a t-distribution gives previous results", {
  skip_on_cran()

  set.seed(26443522)
  t_fit <- suppressWarnings(PF_EM(
    formula = survival::Surv(start, stop, event) ~ ddFixed(group),
    data = head_neck_cancer, id = head_neck_cancer$id, model = "logit",
    by = 1, Q_0 = 1, Q = .1^2,
    control = PF_control(
      N_fw_n_bw = 500, N_smooth = 500, N_first = 1000, n_max = 1,
      Q_tilde = diag(0, 1), nu = 5L),
    max_T = 45, type = "RW"))

  expect_lte(
    t_fit$effective_sample_size$forward_clouds[1],
    length(t_fit$clouds$forward_clouds[[1L]]$weights))
  expect_lte(
    tail(t_fit$effective_sample_size$backward_clouds, 1),
    length(tail(t_fit$clouds$backward_clouds, 1)[[1L]]$weights))

  expect_known_value(t_fit[!names(t_fit) %in% c("clouds", "call", "terms")],
                     file = "pf_t_dist_proposal.RDS")
})

test_that("`get_Q_0` returns a real matrix also when `Fmat` has a complex eigendecomposition", {
  Fmat <- structure(c(0.657881453882877, 0.29770504266233, -0.0894041006717788,
                      0.361585659007044), .Dim = c(2L, 2L))
  Qmat <- structure(c(0.037357449287092, -0.00781178357213716, -0.00781178357213716,
                      0.0205796045976825), .Dim = c(2L, 2L))

  out <- get_Q_0(Qmat = Qmat, Fmat = Fmat)
  expect_true(!is.complex(out))
  expect_equal(
    out,
    structure(c(0.065269831500739, 0.00500931948325941, 0.0050093194832594,
                0.031570480294995), .Dim = c(2L, 2L)))
})

test_that("'PF_forward_filter' gives the same as 'PF_EM' when it should", {
  skip_on_cran()

  #####
  # random walk model
  set.seed(seed <- 56219373)
  pf_em_res <- suppressWarnings(PF_EM(
    survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer, type = "RW", a_0 = c(0, 0),
    fixed_effects = numeric(),
    by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
    control = PF_control(
      N_fw_n_bw = 25, N_smooth = 25, N_first = 25, n_max = 1,
      n_threads = 1), max_T = 30))

  set.seed(seed)
  filter_res <- PF_forward_filter(
    survival::Surv(start, stop, event) ~ group, N_fw = 25, N_first = 25,
    data = head_neck_cancer, type = "RW", a_0 = c(0, 0),
    fixed_effects = numeric(), by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
    control = PF_control(
      N_fw_n_bw = 25, N_smooth = 25, N_first = 25, n_max = 1,
      n_threads = 1), max_T = 30)

  expect_equal(filter_res$forward_clouds, pf_em_res$clouds$forward_clouds)

  pf_em_res[c("a_0", "Q")] <- list(c(0, 0), diag(0.1, 2))
  filter_res <- PF_forward_filter(pf_em_res, N_fw = 25, N_first = 25)
  expect_equal(filter_res$forward_clouds, pf_em_res$clouds$forward_clouds)

  #####
  # dense VAR model with fixed and random arguments
  set.seed(seed)
  pf_em_res <- suppressWarnings(PF_EM(
    fixed = survival::Surv(start, stop, event) ~ 1,
    random = ~ group,
    data = head_neck_cancer, type = "VAR", a_0 = c(0, 0),
    fixed_effects = -2, Fmat = matrix(c(.33, .1, .1, .33), 2),
    by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
    control = PF_control(
      N_fw_n_bw = 25, N_smooth = 25, N_first = 25, n_max = 1,
      n_threads = 1), max_T = 30))

  set.seed(seed)
  filter_res <- PF_forward_filter(
    head_neck_cancer, N_fw = 25, N_first = 25,
    fixed = survival::Surv(start, stop, event) ~ 1,
    random = ~ group,
    type = "VAR", a_0 = c(0, 0),
    fixed_effects = -2, by = 1, Q_0 = diag(1, 2), Q = diag(0.1, 2),
    Fmat = matrix(c(.33, .1, .1, .33), 2),
    control = PF_control(
      N_fw_n_bw = 25, N_smooth = 25, N_first = 25, n_max = 1,
      n_threads = 1), max_T = 30)
  expect_equal(filter_res$forward_clouds, pf_em_res$clouds$forward_clouds)

  pf_em_res[c("fixed_effects", "Q", "F")] <-
    list(-2, diag(0.1, 2), matrix(c(.33, .1, .1, .33), 2))
  filter_res <- PF_forward_filter(pf_em_res, N_fw = 25, N_first = 25)
  expect_equal(filter_res$forward_clouds, pf_em_res$clouds$forward_clouds)

  #####
  # non-dense VAR model
  G <- matrix(c(1, 0, 0, 1), ncol = 1)
  J <- matrix(c(1, 1), ncol = 1)
  K <- matrix(1, 1)

  set.seed(seed)
  pf_em_res <- suppressWarnings(PF_EM(
    survival::Surv(start, stop, event) ~ group + ddFixed_intercept(TRUE),
    data = head_neck_cancer, type = "VAR", a_0 = c(0, 0),
    fixed_effects = -2, G = G, J = J, K = K, psi = log(.1),
    phi = 1, theta = .9, by = 1, Q_0 = diag(1, 2),
    control = PF_control(
      N_fw_n_bw = 25, N_smooth = 25, N_first = 25, n_max = 1,
      n_threads = 1), max_T = 30))

  set.seed(seed)
  filter_res <- PF_forward_filter(
    head_neck_cancer, N_fw = 25, N_first = 25,
    formula = survival::Surv(start, stop, event) ~
      group + ddFixed_intercept(TRUE),
    type = "VAR", a_0 = c(0, 0), G = G, J = J, K = K, phi = 1, psi = log(.1),
    theta = .9, fixed_effects = -2, by = 1, Q_0 = diag(1, 2),
    control = PF_control(
      N_fw_n_bw = 25, N_smooth = 25, N_first = 25, n_max = 1,
      n_threads = 1), max_T = 30)
  expect_equal(filter_res$forward_clouds, pf_em_res$clouds$forward_clouds)

  pf_em_res[c("theta", "phi", "psi", "fixed_effects")] <-
    list(.9, 1, log(.1), -2)
  filter_res <- PF_forward_filter(pf_em_res, N_fw = 25, N_first = 25)
  expect_equal(filter_res$forward_clouds, pf_em_res$clouds$forward_clouds)
})

if(dir.exists("pf-internals"))
  invisible(sapply(list.files("pf-internals", full.names = TRUE), source))

sim_test_pf_internal <- function(n, p, fam, state){
  X <- matrix(runif(n * p, -1, 1), nrow = p)
  offset <- runif(n, -1, 1)

  tstart <- numeric(n)
  tstop <- tstart + 1 + runif(n, max = 3)

  bin_start <- 1
  bin_stop  <- 2

  eta <- drop(state %*% X) + offset
  if(fam == "binomial"){
    y <- 1/(1 + exp(-eta)) > runif(n)

  } else if(fam == "cloglog"){
    y <- -expm1(-exp(eta)) > runif(n)

  } else if(fam == "poisson"){ # TODO: poisson is missleading
    y <- rexp(length(eta), exp(eta))
    dt <- pmin(tstop, bin_stop) - pmax(tstart, bin_start)
    is_event <- y <= dt
    tstop[is_event] <- pmax(tstart[is_event], bin_start) + y[is_event]
    y <- Surv(pmin(y, dt), is_event)

  }

  list(y = y, eta = eta, bin_start = bin_start,
       bin_stop = bin_stop, tstart = tstart, tstop = tstop,
       offset = offset, X = X)
}

test_that("'state_fw' gives correct results", {
  skip_if(!dir.exists("pf-internals"))
  skip_if_not_installed("mvtnorm")

  parent <- c(-.25, -.1)
  parent1 <- c(.3, .18)
  chi <- c(.2, -.05)
  chi1 <- c(-.2, -.07)
  F. <- matrix(c(.9, .4, 0, .8), nrow = 2)
  Q <- matrix(c(.8, .4, .4, .4), nrow = 2)

  obj <- fw(parent = parent, F. = F., Q = Q)

  cpp_out <- check_state_fw(
    parent = parent, parent1 = parent1, child = chi, child1 = chi1,
    F = F., Q = Q)

  require(mvtnorm)
  expect_equal(dmvnorm(chi, F. %*% parent, Q, TRUE), cpp_out$log_dens_func)
  expect_true(cpp_out$is_mvn)
  expect_true(!cpp_out$is_grad_z_hes_const)
  expect_equal(cpp_out$dim, length(parent))

  expect_equal(obj$f(chi), cpp_out$log_dens)
  expect_equal(obj$f(chi1), cpp_out$log_dens1)

  expect_equal(obj$deriv(chi), cpp_out$gradient)
  expect_equal(obj$deriv(chi1), cpp_out$gradient1)

  expect_equal(obj$deriv_z(parent), cpp_out$gradient_zero)
  expect_equal(obj$deriv_z(parent1), cpp_out$gradient_zero1)

  expect_equal(obj$n_hessian(chi), cpp_out$neg_Hessian)
  expect_equal(obj$n_hessian(chi1), cpp_out$neg_Hessian1)
})

test_that("'state_bw' gives correct results", {
  skip_if(!dir.exists("pf-internals"))
  skip_if_not_installed("mvtnorm")

  parent <- c(-.25, -.1)
  parent1 <- c(.3, .18)
  chi <- c(.2, -.05)
  chi1 <- c(-.2, -.07)
  F. <- matrix(c(.9, .4, 0, .8), nrow = 2)
  Q <- matrix(c(.8, .4, .4, .4), nrow = 2)

  obj <- bw(child = chi, F. = F., Q = Q)

  cpp_out <- check_state_bw(
    parent = parent, parent1 = parent1, child = chi, child1 = chi1,
    F = F., Q = Q)

  require(mvtnorm)
  expect_equal(dmvnorm(chi, F. %*% parent, Q, TRUE), cpp_out$log_dens_func)
  expect_true(cpp_out$is_mvn)
  expect_true(!cpp_out$is_grad_z_hes_const)
  expect_equal(cpp_out$dim, length(chi))
  expect_equal(obj$f(parent), cpp_out$log_dens)
  expect_equal(obj$f(parent1), cpp_out$log_dens1)

  expect_equal(obj$deriv(parent), cpp_out$gradient)
  expect_equal(obj$deriv(parent1), cpp_out$gradient1)

  expect_equal(obj$deriv_z(chi), cpp_out$gradient_zero)
  expect_equal(obj$deriv_z(chi1), cpp_out$gradient_zero1)

  expect_equal(obj$n_hessian(chi), cpp_out$neg_Hessian)
  expect_equal(obj$n_hessian(chi1), cpp_out$neg_Hessian1)
})

test_that("'artificial_prior' gives correct results", {
  skip_if(!dir.exists("pf-internals"))
  skip_if_not_installed("mvtnorm")

  state <- c(-.25, -.1)
  F. <- matrix(c(.9, .4, 0, .8), nrow = 2)
  Q <- matrix(c(.8, .4, .4, .4), nrow = 2)
  m_0 <- c(.2, .3)

  prio <- prior(F. = F., Q = Q, Q_0 = Q, mu_0 = m_0)

  t1 <- 1L
  t2 <- 4L
  t3 <- 20L
  ts <- c(t1, t2, t3)

  cpp_out <- check_artificial_prior(
    state = state, F = F., Q = Q, m_0 = m_0, Q_0 = Q, t1 = t1, t2 = t2,
    t3 = t3)

  for(i in seq_along(cpp_out)){
    cpp_o <- cpp_out[[i]]
    p <- prio(ts[i])

    expect_true(cpp_o$is_mvn)
    expect_true(cpp_o$is_grad_z_hes_const)
    expect_equal(cpp_o$dim, length(state))
    expect_equal(cpp_o$log_dens, p$f(state))
    expect_equal(cpp_o$gradient, p$deriv(state))
    expect_equal(cpp_o$gradient_zero, p$deriv_z(state))
    expect_equal(cpp_o$neg_Hessian, p$n_hessian(state))
  }
})

test_that("'observational_cdist' gives correct results", {
  skip_if(!dir.exists("pf-internals"))
  skip_if_not_installed("mvtnorm")

  p <- 3
  state1 <- 1:p - p / 2
  state2 <- state1 + 1
  set.seed(1)

  for(l in list(
    list("binomial", binomial()), list("cloglog", binomial("cloglog")),
    list("poisson", "exponential"))){
    fam <- l[[1]]
    o <- sim_test_pf_internal(n = 50L, p = 3L, fam = fam, state = state1)

    obj <- with(o, odist(y, t(X), offset, l[[2]]))

    y_use <- if(fam == "poisson") o$y[, 2] else o$y
    cpp_out <- with(o, check_observational_cdist(
      X = X, is_event = y_use, offsets = offset, tstart = tstart,
      tstop = tstop, bin_start = bin_start, bin_stop = bin_stop,
      multithreaded = FALSE, fam = fam, state = state1,
      state1 = state2))

    cpp_out_mult <- with(o, check_observational_cdist(
      X = X, is_event = y_use, offsets = offset, tstart = tstart,
      tstop = tstop, bin_start = bin_start, bin_stop = bin_stop,
      multithreaded = TRUE, fam = fam, state = state1,
      state1 = state2))
    expect_equal(cpp_out, cpp_out_mult)

    expect_true(!cpp_out$is_mvn)
    expect_equal(cpp_out$dim, length(state1))

    expect_equal(cpp_out$log_dens, obj$f(state1))
    expect_equal(cpp_out$log_dens1, obj$f(state2))

    expect_equal(drop(cpp_out$gradient), obj$deriv(state1))
    expect_equal(drop(cpp_out$gradient1), obj$deriv(state2))

    expect_equal(cpp_out$neg_Hessian, obj$n_hessian(state1))
    expect_equal(cpp_out$neg_Hessian1, obj$n_hessian(state2))
  }
})

test_that("combining forward and backwards works", {
  skip_if(!dir.exists("pf-internals"))
  skip_if_not_installed("mvtnorm")

  p1 <- c(-.5, 0)
  p2 <- c(-.2, .2)
  c1 <- p2
  c2 <- p1
  x <- c(-.33, .33)
  F. <- matrix(c(.9, .4, 0, .8), nrow = 2)
  Q <- matrix(c(.8, .4, .4, .4), nrow = 2)

  for(nu in c(-1L, 1L, 4L)){
    cpp_out <- check_fw_bw_comb(
      F = F., Q = Q, parent = p1, parent1 = p2,
      grand_child = c1, grand_child1 = c2, x = x, nu = nu)

    #####
    Sig_org <- Sig <- solve(solve(Q) + crossprod(F., solve(Q, F.)))
    m1 <- Sig %*% (solve(Q, F. %*% p1) + crossprod(F., solve(Q, c1)))

    if(nu <= 2L){
      lfunc <- dmvnorm

    } else {
      Sig <- Sig * (nu - 2L) / nu
      lfunc <- function(x, m, s, log)
        dmvt(x = x, delta = m, sigma = s, df = nu, log = TRUE)

    }

    cpp_1 <- cpp_out[[1L]]
    expect_equal(cpp_1$mean, m1)
    expect_equal(cpp_1$covar, Sig)
    require(mvtnorm)
    expect_equal(cpp_1$log_density, lfunc(x, m1, Sig, log = TRUE))

    #####
    m2 <- Sig_org %*% (solve(Q, F. %*% p2) + crossprod(F., solve(Q, c2)))

    cpp_2 <- cpp_out[[2L]]
    expect_equal(cpp_2$mean, m2)
    expect_equal(cpp_2$covar, Sig)
    expect_equal(cpp_2$log_density, lfunc(x, m2, Sig, log = TRUE))
  }
})

test_that("combining prior and backwards works", {
  skip_if(!dir.exists("pf-internals"))
  skip_if_not_installed("mvtnorm")

  p <- c(-.5, 0)
  c1 <- c(-.4, 1)
  c2 <- c(.2, .2)
  m_0 <- c(0, .4)
  F. <- matrix(c(.9, .4, 0, .8), nrow = 2)
  Q <- matrix(c(.8, .4, .4, .4), nrow = 2)
  Q_0 <- Q * 1.1
  t1 <- 1L
  t2 <- 10L

  bwo <- bw(child = c1, F. = F., Q = Q)
  prio <- prior(F. = F., Q = Q, Q_0 = Q_0, mu_0 = m_0)

  cpp_out <- check_prior_bw_comb(
    F = F., Q = Q, m_0 = m_0, Q_0 = Q_0, child = c1, child1 = c2,
    parent = p, t1 = t1, t2 = t2)

  for(t. in c(t2, t1)){
    co <- cpp_out[[as.character(t.)]]
    app <- approximator(prio(t.), bwo, start = c1)

    o1 <- app(NULL, c1)
    o2 <- app(NULL, c2)
    expect_equal(drop(co$mean1), o1$mu)
    expect_equal(drop(co$mean2), o2$mu)

    expect_equal(co$covar1, o1$Sig)
    expect_equal(co$covar2, o2$Sig)

    expect_equal(co$log_dens1, o1$p_ldens(p))
    expect_equal(co$log_dens2, o2$p_ldens(p))
  }
})

test_that("mode approximations give expected result", {
  skip_if(!dir.exists("pf-internals"))
  skip_if_not_installed("mvtnorm")
  skip_if_not_installed("nloptr")

  p <- c(-.5, 0)
  c1 <- c(-.4, 1)
  c2 <- c(.2, .2)
  m_0 <- c(0, .4)
  F. <- matrix(c(.9, .4, 0, .8), nrow = 2)
  Q <- matrix(c(.8, .4, .4, .4), nrow = 2)
  Q_0 <- Q * 1.1
  t1 <- 3L

  bwo <- bw(child = c1, F. = F., Q = Q)
  prio <- prior(F. = F., Q = Q, Q_0 = Q_0, mu_0 = m_0)(t1)

  for(l in list(
    list("poisson", "exponential"),
    list("cloglog", binomial("cloglog")), list("binomial", binomial()))){
    fam <- l[[1]]
    set.seed(1)
    o <- sim_test_pf_internal(n = 50L, p = 2L, fam = fam, state = p)

    y_dist <- with(o, odist(y, t(X), offset, l[[2]]))

    app1 <- approximator(bwo, prio, start = c1)
    app <- approximator(bwo, prio, y_dist, start = app1(c1, NULL)$mu)
    o1 <- app(c1, NULL)
    o2 <- app(c2, NULL)
    Q_xtra <- Q
    Q_xtra[] <- 0

    y_use <- if(fam == "poisson") o$y[, 2] else o$y
    cpp_out <- with(o, check_prior_bw_state_comb(
      X = X, is_event = y_use, offsets = offset, tstart = tstart,
      tstop = tstop, bin_start = bin_start, bin_stop = bin_stop, fam = fam,
      F = F., Q = Q, Q_0 = Q_0, m_0 = m_0, child = c1, child1 = c2, parent = p,
      t1 = t1, Q_xtra = Q_xtra))

    expect_equal(drop(cpp_out$mean1), o1$mu)
    expect_equal(drop(cpp_out$mean2), o2$mu)

    expect_equal(cpp_out$covar1, o1$Sig)
    expect_equal(cpp_out$covar2, o2$Sig)

    expect_equal(cpp_out$log_dens1, o1$p_ldens(p))
    expect_equal(cpp_out$log_dens2, o2$p_ldens(p))
  }

  #####
  # with extra covariance matrix term
  Q_xtra <- diag(sqrt(.3), 2)
  set.seed(1)
  o <- sim_test_pf_internal(n = 50L, p = 2L, fam = fam, state = p)

  y_dist <- with(o, odist(y, t(X), offset, l[[2]]))

  app1 <- approximator(bwo, prio, start = c1)
  app <- approximator(bwo, prio, y_dist, start = app1(c1, NULL)$mu)
  o1 <- app(c1, NULL)
  o2 <- app(c2, NULL)

  y_use <- if(fam == "poisson") o$y[, 2] else o$y
  cpp_out <- with(o, check_prior_bw_state_comb(
    X = X, is_event = y_use, offsets = offset, tstart = tstart,
    tstop = tstop, bin_start = bin_start, bin_stop = bin_stop, fam = fam,
    F = F., Q = Q, Q_0 = Q_0, m_0 = m_0, child = c1, child1 = c2, parent = p,
    t1 = t1, Q_xtra = Q_xtra))

  expect_equal(drop(cpp_out$mean1), o1$mu)
  expect_equal(drop(cpp_out$mean2), o2$mu)

  Q_res <- o1$Sig + Q_xtra
  expect_equal(cpp_out$covar1, Q_res)
  expect_equal(cpp_out$covar2, Q_res)

  expect_equal(cpp_out$log_dens1, dmvnorm(p, o1$mu, Q_res, TRUE))
  expect_equal(cpp_out$log_dens2, dmvnorm(p, o2$mu, Q_res, TRUE))

  #####
  # with extra covariance matrix term and t-distribution
  nu <- 4
  cpp_out <- with(o, check_prior_bw_state_comb(
    X = X, is_event = y_use, offsets = offset, tstart = tstart,
    tstop = tstop, bin_start = bin_start, bin_stop = bin_stop, fam = fam,
    F = F., Q = Q, Q_0 = Q_0, m_0 = m_0, child = c1, child1 = c2, parent = p,
    t1 = t1, Q_xtra = Q_xtra, nu = nu))

  expect_equal(drop(cpp_out$mean1), o1$mu)
  expect_equal(drop(cpp_out$mean2), o2$mu)

  Q_res <- (o1$Sig + Q_xtra) * (nu - 2) / nu
  expect_equal(cpp_out$covar1, Q_res)
  expect_equal(cpp_out$covar2, Q_res)

  expect_equal(cpp_out$log_dens1, dmvt(p, o1$mu, Q_res, nu, TRUE))
  expect_equal(cpp_out$log_dens2, dmvt(p, o2$mu, Q_res, nu, TRUE))

  #####
  # with extra covariance matrix term, scaled and t-distribution
  nu <- 4
  covar_fac <- 1.5
  cpp_out <- with(o, check_prior_bw_state_comb(
    X = X, is_event = y_use, offsets = offset, tstart = tstart,
    tstop = tstop, bin_start = bin_start, bin_stop = bin_stop, fam = fam,
    F = F., Q = Q, Q_0 = Q_0, m_0 = m_0, child = c1, child1 = c2, parent = p,
    t1 = t1, Q_xtra = Q_xtra, nu = nu, covar_fac = covar_fac))

  expect_equal(drop(cpp_out$mean1), o1$mu)
  expect_equal(drop(cpp_out$mean2), o2$mu)

  Q_res <- (o1$Sig + Q_xtra) * (nu - 2) / nu * covar_fac
  expect_equal(cpp_out$covar1, Q_res)
  expect_equal(cpp_out$covar2, Q_res)

  expect_equal(cpp_out$log_dens1, dmvt(p, o1$mu, Q_res, nu, TRUE))
  expect_equal(cpp_out$log_dens2, dmvt(p, o2$mu, Q_res, nu, TRUE))
})

test_that("'PF_get_score_n_hess' gives the same as an R implementation", {
  skip_on_cran()

  .lung <- lung[!is.na(lung$ph.ecog), ]
  # standardize
  .lung$age <- scale(.lung$age)

  # fit
  set.seed(43588155)
  pf_fit <- suppressWarnings(PF_EM(
    Surv(time, status == 2) ~ ddFixed(ph.ecog) + ddFixed(age) + age,
    data = .lung, by = 50, id = 1:nrow(.lung), Q_0 = diag(1, 2),
    Fmat = structure(c(0.96, -0.00045, 0.066, 0.37), .Dim = c(2L, 2L)),
    Q = structure(c(0.071, -0.055, -0.055, 0.058), .Dim = c(2L, 2L)),
    max_T = 800, fixed_effects = c(0.36, 0.1), type = "VAR",
    control = PF_control(
      N_fw_n_bw = 500, N_first = 2500, N_smooth = 100,
      n_max = 1, eps = .001, Q_tilde = diag(.2^2, 2), est_a_0 = FALSE,
      n_threads = max(parallel::detectCores(logical = FALSE), 1))))

  test_logit_func <- function(fit){
    org_cl <- fit$call
    ma <- match(
      c("formula", "data", "by", "max_T", "id", "trace", "model", "order",
        "fixed", "fixed", "random"), names(org_cl), nomatch = 0L)
    sta_arg_call <- org_cl[c(1L, ma)]

    # add defaults where needed
    def_args <- c("trace", "model", "fixed", "random", "order")
    if(any(is_missing <- !def_args %in% names(sta_arg_call)))
      sta_arg_call[def_args[is_missing]] <- formals(PF_EM)[def_args[is_missing]]
    sta_arg_call[[1L]] <- quote(dynamichazard:::.get_PF_static_args)

    static_args <- eval(sta_arg_call)
    fw <- fit$clouds$forward_clouds

    scores_old <- hess_old <- NULL
    k <- length(fit$F)
    K <- solve(fit$Q)
    K1_2 <- K / 2
    for(i in seq_len(length(fw) - 1L) + 1L){
      fw_i <- fw[[i]]
      ps <- fw_i$states
      par_idx <- fw_i$parent_idx
      parents <- fw[[i - 1L]]$states[, par_idx]

      R_i <- static_args$risk_obj$risk_sets[[i - 1L]]
      X_i <- static_args$X[, R_i, drop = FALSE]
      Z_i <- static_args$fixed_terms[, R_i, drop = FALSE]
      y_i <- static_args$risk_obj$is_event_in[R_i] == (i - 2L)

      eta <- drop(crossprod(Z_i, fit$fixed_effects))

      nf <- NROW(Z_i)
      score_dim <- nf + k * 2L
      scores <- matrix(0., ncol(ps), score_dim)
      hess <- array(0, c(score_dim, score_dim, ncol(ps)))

      #####
      # state equation
      r0 <- ps - fit$F %*% parents
      r1 <- solve(fit$Q, r0)
      for(j in 1:ncol(ps)){
        scores[j, nf + 1:k] <- tcrossprod(parents[, j], r1[, j])
        scores[j, nf + 1:k + k] <-
          solve(fit$Q, tcrossprod(r0[, j] * .5, r1[, j]) - diag(.5, k / 2L))

        hess[nf + 1:k, nf + 1:k, j] <- - kronecker(K, tcrossprod(parents[, j]))
        off_block <- - kronecker(K, tcrossprod(r1[, j], parents[, j]))
        hess[nf + 1:k + k, nf + 1:k    , j] <- off_block
        hess[nf + 1:k    , nf + 1:k + k, j] <- t(off_block)
        hess[nf + 1:k + k, nf + 1:k + k, j] <-
          - kronecker(K, tcrossprod(r1[, j]) - K1_2)
      }

      #####
      # observational equation
      a_obs <- apply(ps, 2, function(x){
        eta <- eta + drop(crossprod(X_i, x))
        exp_eta <- exp(eta)
        drop(crossprod((exp_eta * (y_i - 1) + y_i) / (1 + exp_eta), t(Z_i)))
      })
      scores[, 1:nf] <- if(!is.matrix(a_obs)) as.matrix(a_obs) else t(a_obs)
      B_obs <- apply(ps, 2, function(x){
        eta <- eta + drop(crossprod(X_i, x))
        exp_eta <- exp(eta)
        crossprod(-t(Z_i) * exp(eta) / (1 + exp(eta))^2, t(Z_i))
      })
      hess[1:nf, 1:nf, ] <- array(B_obs, dim = c(nf, nf, dim(hess)[[3]]))

      if(!is.null(scores_old))
        scores <- scores + scores_old[par_idx, ]
      if(!is.null(hess_old))
        hess <- hess + hess_old[, , par_idx]

      scores_old <- scores
      hess_old   <- hess
    }

    ws <- drop(tail(fw, 1)[[1L]]$weights)
    score <- colSums(scores * ws)

    info_obj <- matrix(0., NCOL(scores), NCOL(scores))
    for(i in seq_along(ws))
      info_obj <-
        info_obj + ws[i] * (tcrossprod(scores[i, ]) + hess[, , i])

    dfix <- NROW(Z_i)
    n_rng <- NROW(X_i)
    out_dim <- dfix + n_rng * n_rng + (n_rng * (n_rng + 1L)) / 2L
    trans_mat <- matrix(0, out_dim, length(score))

    ifix <- 1:dfix
    trans_mat[ifix, ifix]          <- diag(dfix)
    idF  <- dfix + 1:(n_rng * n_rng)
    trans_mat[idF, idF]            <- dynamichazard:::.get_cum_mat(n_rng, n_rng)
    idQ <- dfix + n_rng * n_rng + 1:((n_rng * (n_rng + 1L)) / 2L)
    trans_mat[idQ, -c(ifix, idF)] <- t(dynamichazard:::.get_dup_mat(n_rng))

    obs_info = tcrossprod(score) - info_obj
    list(
      score = drop(trans_mat %*% score),
      obs_info = tcrossprod(trans_mat %*% obs_info, trans_mat))
  }

  test_out <- test_logit_func(pf_fit)
  expect_output(func_out <- PF_get_score_n_hess(pf_fit),
                regexp = "Using '.lung' as the 'data' argument", fixed = TRUE)

  expect_equal(test_out, func_out$get_get_score_n_hess(only_score = FALSE),
               check.attributes = FALSE)

  # only score
  test_out$obs_info[] <- NA_real_
  expect_equal(test_out, func_out$get_get_score_n_hess(only_score = TRUE),
               check.attributes = FALSE)
})

test_that("warns with randow walk model when covariates are in both fixed and random", {
  skip_on_cran()

  set.seed(43588155)
  expect_warning(
    pf_fit <- PF_EM(
      fixed = Surv(tstart, tstop, death == 2) ~
        age + edema  + log(albumin) + log(protime) + log(bili),
      random = ~ log(bili),
      data = pbc2, Q_0 = diag(1, 2), Q = diag(1, 2),
      by = 50L, id = pbc2$id, a_0 = c(0, 0),
      fixed_effects = c(-15, 0, 0, 0, 0, 0),
      model = "exponential", max_T = 100L,
      control = PF_control(
        N_fw_n_bw = 50, N_smooth = 50, N_first = 50, eps = 1e-3,
        method = "AUX_normal_approx_w_cloud_mean", est_a_0 = FALSE,
        n_max = 1L)),
    regexp = "The following terms are in both 'fixed' and 'random' with 'type = \"RW\"':  'intercept', 'log(bili)'",
    fixed = TRUE)
})

test_that("works when there are time periods outside where we have data and with a single fixed and random effect", {
  skip_on_cran()
  dir <- "local_tests"
  skip_if_not(dir.exists(file.path("previous_results", dir)))

  cap_ti <- 20L
  head_neck_cancer <- subset(head_neck_cancer, stop < cap_ti)

  set.seed(43588155)
  suppressWarnings(pf_fit <- PF_EM(
    fixed = Surv(start, stop, event) ~ 1, random = ~ I((group == 2) + 0.) - 1,
    data = head_neck_cancer, Q_0 = diag(1, 1), Q = diag(0.0923169242234624, 1),
    fixed_effects = -2.68255823834327,
    a_0 = -0.745106571679811,
    by = 1L, id = head_neck_cancer$id, type = "RW",
    model = "cloglog", max_T = cap_ti + 6L,
    control = PF_control(
      N_fw_n_bw = 400, N_smooth = 1000L, N_first = 1000L,
      method = "AUX_normal_approx_w_cloud_mean", est_a_0 = FALSE,
      n_max = 1L)))

  expect_known_value(
    get_cloud_means(pf_fit), file.path(dir, "pf-rw-w-no-obs-periods.RDS"))
})

test_that("'fix_seed = FALSE' yields different results at each call", {
  skip_on_cran()

  set.seed(seed <- 1L)
  same_seed_call <- suppressWarnings(
    PF_EM(fixed = Surv(stop, event) ~ group, random = ~ 1, model = "logit",
          max_T = 30L, by = 1L, data = head_neck_cancer,
          Q = as.matrix(1), Fmat = as.matrix(.9),
          type = "VAR", control = PF_control(
            N_fw_n_bw = 100L, N_smooth = 100L, N_first = 100L, n_max = 3L,
            fix_seed = TRUE)))

  call_2 <- same_seed_call$call
  call_2[["control"]]$fix_seed <- FALSE
  set.seed(seed)
  diff_seed <- suppressWarnings(eval(call_2, environment()))

  # just check one element
  expect_equal(diff_seed$EM_ests     $fixed_effects[1, ],
               same_seed_call$EM_ests$fixed_effects[1, ])

  expect_true(!isTRUE(all.equal(
    diff_seed$EM_ests     $fixed_effects[-1, ],
    same_seed_call$EM_ests$fixed_effects[-1, ])))

  # check particle filter
  f1 <- PF_forward_filter(diff_seed, N_fw = 100L, N_first = 100L)
  f2 <- PF_forward_filter(diff_seed, N_fw = 100L, N_first = 100L)
  expect_true(!isTRUE(all.equal(f1, f2)))

  # check score and Hessian method
  o <- PF_get_score_n_hess(diff_seed)
  o$set_n_particles(N_fw = 100L, N_first = 100L)
  expect_true(!isTRUE(all.equal(o$run_particle_filter(),
                                o$run_particle_filter())))

  o <- PF_get_score_n_hess(diff_seed, use_O_n_sq = TRUE)
  o$set_n_particles(N_fw = 50L, N_first = 50L)
  expect_true(!isTRUE(all.equal(o$get_get_score_n_hess(),
                                o$get_get_score_n_hess())))
})
