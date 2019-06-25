context("Testing examples in man files")

test_that("ddhazard help page examples gives the same results", {
  pbc <- pbc_org
  fit <- ddhazard(
    Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
    Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
    control = ddhazard_control(method = "GMA"))

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
   control = ddhazard_control(method = "GMA"),
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
    control = ddhazard_control(method = "GMA"))

  resids <- residuals(fit, type = "pearson")$residuals
  expect_known_value(resids[1:2], "ddhazard_help_page_residuals.RDS")
})

test_that("hatvalues.ddhazard help page examples gives the same results", {
  pbc <- pbc_org
  fit <- ddhazard(
    Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3000,
    Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 100,
    control = ddhazard_control(method = "GMA"))
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
    control = ddhazard_control(method = "GMA"))
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
get_means <- function(result)
  lapply(
    result[c("forward_clouds", "backward_clouds", "smoothed_clouds")],
    function(clouds)
      do.call(rbind, sapply(clouds, function(row){
        colSums(t(row$states) * drop(row$weights))
      }, simplify = FALSE)))

test_that("PF_EM help page example runs and gives previous computed results", {
  skip_on_cran()

  .lung <- lung[!is.na(lung$ph.ecog), ]
  # standardize
  .lung$age <- scale(.lung$age)

  # fit
  set.seed(43588155)
  pf_fit <- suppressWarnings(PF_EM(
    Surv(time, status == 2) ~ ddFixed(ph.ecog) + age,
    data = .lung, by = 50, id = 1:nrow(.lung),
    Q_0 = diag(1, 2), Q = diag(.5^2, 2),
    max_T = 800,
    control = PF_control(
      N_fw_n_bw = 500, N_first = 2500, N_smooth = 5000,

      # take less its than in the man page
      n_max = 5,

      eps = .001, Q_tilde = diag(.2^2, 2), est_a_0 = FALSE,
      n_threads = max(parallel::detectCores(logical = FALSE), 1))))

  if(dir.exists("previous_results/local_tests"))
    # tmp <- readRDS("previous_results/local_tests/survival_lung_example.RDS")
    expect_known_value(
      pf_fit[!names(pf_fit) %in% c("call", "control")],
      "local_tests/survival_lung_example.RDS", tolerance = 1.49e-08)
  expect_equal(
    pf_fit$control$n_threads, max(parallel::detectCores(logical = FALSE), 1))

  # tmp <- readRDS("previous_results/survival_lung_example_cloud_means.RDS")
  expect_known_value(
    get_means(pf_fit$clouds), "survival_lung_example_cloud_means.RDS",
    tolerance = 1.49e-08)

  expect_no_error(plot(pf_fit, cov_index = 1))
  expect_no_error(plot(pf_fit, cov_index = 2))
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
    control = ddhazard_control(eps = 1e-5, NR_eps = 1e-4, n_max = 1e4))

  expect_known_value(ddfit[!names(ddfit) %in% c("call", "control")],
                     "local_tests/pf_man_2nd_ddfit.RDS")

  # fit model with particle filter
  set.seed(88235076)
  pf_fit <- suppressWarnings(PF_EM(
    Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() + ddFixed(age) +
      ddFixed(edema) + ddFixed(log_albumin) + ddFixed(log_protime) + log_bili,
    pbc2, Q_0 = 2^2, Q = ddfit$Q * 100, # use estimate from before
    by = 100, id = pbc2$id,
    model = "exponential", max_T = 3600,
    control = PF_control(
      N_fw_n_bw = 500, N_smooth = 2500, N_first = 1000, eps = 1e-3,
      method = "AUX_normal_approx_w_cloud_mean", est_a_0 = FALSE,
      Q_tilde = as.matrix(.1^2),
      n_max = 5, # less than in man page
      n_threads = max(parallel::detectCores(logical = FALSE), 1))))

  # tmp <- readRDS("previous_results/local_tests/pf_man_2nd_ppfit.RDS")
  expect_known_value(pf_fit[!names(pf_fit) %in% c("clouds", "call", "control")],
                     "local_tests/pf_man_2nd_ppfit.RDS")
  expect_equal(
    pf_fit$control$n_threads, max(parallel::detectCores(logical = FALSE), 1))
})

test_that("example in 'PF_EM' with gives previous results w/ a few iterations", {
  skip_on_cran()
  skip_if(!dir.exists("previous_results/local_tests"))

  # g groups with k individuals in each
  g <- 3L
  k <- 400L

  # matrices for state equation
  p <- g + 1L
  G <- matrix(0., p^2, 2L)
  for(i in 1:p)
    G[i + (i - 1L) * p, 1L + (i == p)] <- 1L

  theta <- c(.9, .8)
  F. <- matrix(as.vector(G %*% theta), 4L, 4L)

  J <- matrix(0., ncol = 2L, nrow = p)
  J[-p, 1L] <- J[p, 2L] <- 1
  psi <- c(log(c(.3, .1)))

  K <- matrix(0., p * (p - 1L) / 2L, 2L)
  j <- 0L
  for(i in (p - 1L):1L){
    j <- j + i
    K[j, 2L] <- 1
  }
  K[K[, 2L] < 1, 1L] <- 1
  phi <- log(-(c(.8, .3) + 1) / (c(.8, .3) - 1))

  V <- diag(exp(drop(J %*% psi)))
  C <- diag(1, ncol(V))
  C[lower.tri(C)] <- 2/(1 + exp(-drop(K %*% phi))) - 1
  C[upper.tri(C)] <- t(C)[upper.tri(C)]
  Q <- V %*% C %*% V     # covariance matrix in state transition

  Q_0 <- get_Q_0(Q, F.)    # invariant covariance matrix
  beta <- c(rep(-6, g), 0) # all groups have the same long run mean intercept

  # simulate state variables
  set.seed(56219373)
  n_periods <- 300L
  alphas <- matrix(nrow = n_periods + 1L, ncol = p)
  alphas[1L, ] <- rnorm(p) %*% chol(Q_0)
  for(i in 1:n_periods + 1L)
    alphas[i, ] <- F. %*% alphas[i - 1L, ] + drop(rnorm(p) %*% chol(Q))

  alphas <- t(t(alphas) + beta)

  # simulate individuals' outcome
  n_obs <- g * k
  df <- lapply(1:n_obs, function(i){
    # find the group
    grp <- (i - 1L) %/% (n_obs / g) + 1L

    # left-censoring
    tstart <- max(0L, sample.int((n_periods - 1L) * 2L, 1) - n_periods + 1L)

    # covariates
    x <- c(1, rnorm(1))

    # outcome (stop time and event indicator)
    osa <- NULL
    oso <- NULL
    osx <- NULL
    y <- FALSE
    for(tstop in (tstart + 1L):n_periods){
      sigmoid <- 1 / (1 + exp(- drop(x %*% alphas[tstop + 1L, c(grp, p)])))
      if(sigmoid > runif(1)){
        y <- TRUE
        break
      }
      if(.01 > runif(1L) && tstop < n_periods){
        # sample new covariate
        osa <- c(osa, tstart)
        tstart <- tstop
        oso <- c(oso, tstop)
        osx <- c(osx, x[2])
        x[2] <- rnorm(1)
      }
    }

    cbind(
      tstart = c(osa, tstart), tstop = c(oso, tstop),
      x = c(osx, x[2]), y = c(rep(FALSE, length(osa)), y), grp = grp,
      id = i)
  })
  df <- data.frame(do.call(rbind, df))
  df$grp <- factor(df$grp)

  # fit model. Start with "cheap" iterations
  fit <- suppressWarnings(PF_EM(
    fixed = Surv(tstart, tstop, y) ~ x, random = ~ grp + x - 1,
    data = df, model = "logit", by = 1L, max_T = max(df$tstop),
    Q_0 = diag(1.5^2, p), id = df$id, type = "VAR",
    G = G, theta = c(0.905, .292), J = J, psi = log(c(0.261, 0.116)),
    K = K, phi = log(-(c(0.63, -0.015) + 1) / (c(0.63, -0.015) - 1)),
    fixed_effects = c(-5.8, 0.049),
    control = PF_control(
      N_fw_n_bw = 100L, N_smooth = 100L, N_first = 500L,
      method = "AUX_normal_approx_w_cloud_mean",
      nu = 5L, # sample from multivariate t-distribution
      n_max = 5L,  averaging_start = 1L,
      smoother = "Fearnhead_O_N", eps = 1e-4, covar_fac = 1.2,
      n_threads = 4L # depends on your cpu(s)
    )))

  expect_known_value(fit[!names(fit) %in% c("clouds", "call")],
                     "local_tests/pf_man_restrict_fit_ex.RDS")

  # take more iterations with more particles
  cl <- fit$call
  ctrl <- cl[["control"]]
  ctrl[c("N_fw_n_bw", "N_smooth", "N_first", "n_max",
         "averaging_start")] <- list(200L, 1000L, 5000L, 2L, 1L)
  cl[["control"]] <- ctrl
  cl[c("phi", "psi", "theta")] <- list(fit$phi, fit$psi, fit$theta)
  fit_extra <- suppressWarnings(eval(cl))

  expect_known_value(fit_extra[!names(fit_extra) %in% c("clouds", "call")],
                     "local_tests/pf_man_restrict_fit_ex_more.RDS")

  # plot predicted state variables. We do this just to check that it does not
  # fail
  for(i in 1:p){
    plot(fit_extra, cov_index = i)
    abline(h = 0, lty = 2)
    lines(1:nrow(alphas) - 1, alphas[, i] - beta[i], lty = 3)
  }
})

test_that("`PF_forward_filter` the results stated in the comments and does not alter the .Random.seed as stated on the help page", {
  skip_on_cran()

  # head-and-neck cancer study data
  # see Efron, B. (1988) doi:10.2307/2288857
  is_censored <- c(
    6, 27, 34, 36, 42, 46, 48:51, 51 + c(15, 30:28, 33, 35:37, 39, 40, 42:45))
  head_neck_cancer <- data.frame(
    id = 1:96,
    stop = c(
      1, 2, 2, rep(3, 6), 4, 4, rep(5, 8),
      rep(6, 7), 7, 8, 8, 8, 9, 9, 10, 10, 10, 11, 14, 14, 14, 15, 18, 18, 20,
      20, 37, 37, 38, 41, 45, 47, 47,
      2, 2, 3, rep(4, 4), rep(5, 5), rep(6, 5),
      7, 7, 7, 9, 10, 11, 12, 15, 16, 18, 18, 18, 21,
      21, 24, 25, 27, 36, 41, 44, 52, 54, 59, 59, 63, 67, 71, 76),
    event = !(1:96 %in% is_censored),
    group = factor(c(rep(1, 45 + 6), rep(2, 45))))

  # fit model
  set.seed(61364778)
  ctrl <- PF_control(
    N_fw_n_bw = 500, N_smooth = 2500, N_first = 2000,
    n_max = 1, # set to one as an example
    n_threads = max(parallel::detectCores(logical = FALSE), 1),
    eps = .001, Q_tilde = as.matrix(.3^2), est_a_0 = FALSE)
  pf_fit <- suppressWarnings(
    PF_EM(
      survival::Surv(stop, event) ~ ddFixed(group),
      data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2, control = ctrl,
      max_T = 30))

  # the log-likelihood in the final iteration
  # dput(logLik(pf_fit))
  expect_equal(
    end_log_like <- logLik(pf_fit), structure(
      -258.731786156198, "P(y_t|y_{1:(t-1)})" = c(
        -6.70941384660467,
        -17.4991371234673, -22.7486543765548, -23.1567710570606, -40.7268735916535,
        -32.547747619028, -12.6814755271764, -10.6604404283416, -10.9186968404081,
        -10.6576119071216, -5.64937811314284, -5.45910379291176, -2.10659339971792,
        -9.7217688765877, -7.55621994298938, -4.86522902610044, -1.65390885383481,
        -7.50079904558733, -1.24740698294264, -6.5341593177666, -4.36317037539044,
        -0.95903332721229, -0.902069557998193, -4.24515671289682, -0.840276896118074,
        -0.678457059772247, -4.20700482272061, -0.614170469993203, -0.69285146373326,
        -0.628205801364764),
      df = 2, nobs = NA_integer_, class = "logLik"))

  # gives the same
  seed_now <- .GlobalEnv$.Random.seed
  fw_ps <- PF_forward_filter(
    survival::Surv(stop, event) ~ ddFixed(group), N_fw = 500, N_first = 2000,
    data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2,
    a_0 = pf_fit$a_0, fixed_effects = -0.5370051,
    control = ctrl, max_T = 30, seed = pf_fit$seed)
  expect_true(isTRUE(all.equal(c(end_log_like), c(logLik(fw_ps)))))
  expect_equal(seed_now, .GlobalEnv$.Random.seed)

  # will differ since we use different number of particles
  fw_ps <- PF_forward_filter(
    survival::Surv(stop, event) ~ ddFixed(group), N_fw = 1000, N_first = 3000,
    data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2,
    a_0 = pf_fit$a_0, fixed_effects = -0.5370051,
    control = ctrl, max_T = 30, seed = pf_fit$seed)
  expect_false(isTRUE(all.equal(c(end_log_like), c(logLik(fw_ps)))))
  expect_equal(seed_now, .GlobalEnv$.Random.seed)

  # will differ since we use the final estimates
  fw_ps <- PF_forward_filter(pf_fit, N_fw = 500, N_first = 2000)
  expect_false(isTRUE(all.equal(c(end_log_like), c(logLik(fw_ps)))))
  expect_equal(seed_now, .GlobalEnv$.Random.seed)

  # should give the same when called again
  fw_ps_2 <- PF_forward_filter(pf_fit, N_fw = 500, N_first = 2000)
  expect_equal(fw_ps, fw_ps_2)

  # but not when we change the seed argument
  runif(100)
  fw_ps_3 <- PF_forward_filter(
    pf_fit, N_fw = 500, N_first = 2000, seed = .Random.seed)
  expect_false(isTRUE(all.equal(fw_ps, fw_ps_3)))
})

test_that("ddsurvcurve manual page examples give the same", {
  pbc <- pbc_org
  temp <- subset(pbc, id <= 312, select=c(id:sex, stage))
  pbc2 <- tmerge(temp, temp, id=id, death = event(time, status))
  pbc2 <- tmerge(pbc2, pbcseq, id = id, bili = tdc(day, bili))

  f1 <- suppressWarnings(ddhazard(
    Surv(tstart, tstop, death == 2) ~ ddFixed(log(bili)), pbc2, id = pbc2$id,
    max_T = 3600, Q_0 = 1, Q = 1e-2, by = 100, model = "exponential",
    control = ddhazard_control(method = "EKF", eps = 1e-1, n_max = 1,
                               fixed_terms_method = "M_step")))

  ddcurve <- ddsurvcurve(f1)
  z <- plot(ddcurve, col = "DarkBlue", lwd = 2)
  expect_known_value(ddcurve, file = "ddsurvcurve-fix-cont.RDS")
  expect_equal(z$S(c((0:3) * 1000)), c(
    1, 0.967806836427118, 0.937426234451691, 0.904842879418892))

  nw <- data.frame(bili = 3)
  z <- lines(ddsurvcurve(f1, new_data = nw), col = "DarkBlue")
  expect_equal(z$S(c((0:3) * 1000)), c(
    1, 0.866618493229879, 0.75375407181175, 0.64567672043475))

  f3 <- suppressWarnings(ddhazard(
    Surv(tstart, tstop, death == 2) ~ log(bili), pbc2, id = pbc2$id,
    max_T = 3600, Q_0 = diag(1, 2), Q = diag(1e-2, 2), by = 100, model = "exponential",
    control = ddhazard_control(method = "EKF", eps = 1e-1, n_max = 1)))

  nw <- data.frame(
    bili = c(2.1, 1.9, 3.3, 3.9, 3.8, 3.6, 4, 4.9, 4.2, 5.7, 10.2),
    tstart = c(0L, 225L, 407L, 750L, 1122L, 1479L, 1849L, 2193L, 2564L, 2913L,
               3284L),
    tstop = c(225L, 407L, 750L, 1122L, 1479L, 1849L, 2193L, 2564L, 2913L,
              3284L, 3600L))
  ddcurve <- ddsurvcurve(f3, new_data = nw, tstart = "tstart", tstop = "tstop")
  z <- lines(ddcurve, "darkorange", lwd = 2)
  expect_equal(z$S(c((0:3) * 1000)), c(
    1, 0.831340175284127, 0.694231628800056, 0.423403100827788))

  ddcurve <- ddsurvcurve(f3, new_data = nw[-(1:5), ], tstart = "tstart",
                         tstop = "tstop")
  z <- lines(ddcurve, lty = 2, lwd = 2)
  expect_equal(z$S(c((2:3) * 1000)), c(0.907583926364646, 0.553523972032861))

  #####
  # example with discrete time model
  h1 <- suppressWarnings(ddhazard(
    Surv(stop, event) ~ group, head_neck_cancer, by = 1, max_T = 45,
    Q_0 = diag(2^2, 2), Q = diag(.01^2, 2), control = ddhazard_control(
      method = "GMA", eps = 1e-1, n_max = 1)))

  nw <- data.frame(group = factor(1, levels = 1:2), tstart = 0, tstop = 30)
  ddcurve <- ddsurvcurve(h1, new_data = nw, tstart = "tstart",
                         tstop = "tstop")
  z <- plot(ddcurve, col = "Darkblue")
  expect_known_value(z, file = "ddsurvcurve-fix-disc-1.RDS")

  nw$group <- factor(2, levels = 1:2)
  ddcurve <- ddsurvcurve(h1, new_data = nw, tstart = "tstart",
                         tstop = "tstop")
  z <- lines(ddcurve, col = "DarkOrange")
  expect_known_value(z, file = "ddsurvcurve-fix-disc-2.RDS")
})

test_that("`get_Q_0` example returns the correct covariance matrix", {
  Fmat <- matrix(c(.8, .4, .1, .5), 2, 2)
  Qmat <- matrix(c( 1, .5, .5,  2), 2)

  x1 <- get_Q_0(Qmat = Qmat, Fmat = Fmat)
  x2 <- Qmat
  for(i in 1:101)
    x2 <- tcrossprod(Fmat %*% x2, Fmat) + Qmat
  expect_equal(x1, x2)
})

test_that("'PF_get_score_n_hess' gives previous results", {
  skip_on_cran()
  #####
  # same as manual page
  library(dynamichazard)
  .lung <- lung[!is.na(lung$ph.ecog), ]
  # standardize
  .lung$age <- scale(.lung$age)

  set.seed(43588155)
  suppressWarnings(pf_fit <- PF_EM(
    fixed = Surv(time, status == 2) ~ ph.ecog + age,
    random = ~ 1, model = "exponential",
    data = .lung, by = 50, id = 1:nrow(.lung),
    Q_0 = as.matrix(1), Q = as.matrix(0.040353), type = "VAR",
    max_T = 800, Fmat = as.matrix(0.77156),
    fixed_effects =
      c(`(Intercept)` = -6.343, ph.ecog = 0.41836, age = 0.10481),
    control = PF_control(
      N_fw_n_bw = 250, N_first = 2000, N_smooth = 500, covar_fac = 1.1,
      nu = 6, n_max = 1L, eps = 1e-4, averaging_start = 200L,
      n_threads = max(parallel::detectCores(logical = FALSE), 1))))

  expect_output(comp_obj <- PF_get_score_n_hess(pf_fit),
                "Using '.lung' as the 'data' argument", fixed = TRUE)
  comp_obj$set_n_particles(N_fw = 1000L, N_first = 1000L)
  comp_obj$run_particle_filter()
  o1 <- comp_obj$get_get_score_n_hess()
  expect_known_value(o1, "PF_get_score_n_hess-help-res.RDS",
                     tolerance = .Machine$double.eps^(1/3))

  # check that we get the same when we only request the score
  o1_score <- comp_obj$get_get_score_n_hess(only_score = TRUE)
  expect_equal(o1$score, o1_score$score)

  expect_output(comp_obj <- PF_get_score_n_hess(pf_fit, use_O_n_sq = TRUE),
                "Using '.lung' as the 'data' argument", fixed = TRUE)
  comp_obj$set_n_particles(N_fw = 250L, N_first = 250L)
  o2 <- comp_obj$get_get_score_n_hess()
  expect_known_value(o2, "PF_get_score_n_hess-help-res-O_N_sq.RDS")

  # check that we get the same when we only request the score
  o2_score <- comp_obj$get_get_score_n_hess(only_score = TRUE)
  expect_equal(o2$score, o2_score$score)

  #####
  # we also run the old 2D test
  set.seed(43588155)
  pf_fit <- suppressWarnings(PF_EM(
    fixed = Surv(time, status == 2) ~ ph.ecog + age,
    random = ~ age, model = "exponential",
    data = .lung, by = 50, id = 1:nrow(.lung), type = "VAR", Q_0 = diag(1, 2),
    fixed_effects = c(-6.3, 0.42, 0.22),
    Q = structure(c(0.068, -0.044, -0.044, 0.039), .Dim = c(2L, 2L)),
    Fmat = structure(c(0.67, 0.23, -0.24, 0.63), .Dim = c(2L, 2L)),
    max_T = 800,
    control = PF_control(
      N_fw_n_bw = 250, N_first = 1000, N_smooth = 500, covar_fac = 1.1,
      nu = 6, n_max = 5L, eps = 1e-5, est_a_0 = FALSE, averaging_start = 100L,
      n_threads = 2L)))

  expect_output(comp_obj <- PF_get_score_n_hess(pf_fit),
                "Using '.lung' as the 'data' argument", fixed = TRUE)
  comp_obj$set_n_particles(N_fw = 1000L, N_first = 1000L)
  comp_obj$run_particle_filter()
  o <- comp_obj$get_get_score_n_hess()
  expect_known_value(o, "PF_get_score_n_hess-help-res-old.RDS")

  # check that we get the same when we only request the score
  o_score <- comp_obj$get_get_score_n_hess(only_score = TRUE)
  expect_equal(o$score, o_score$score)

  expect_output(comp_obj <- PF_get_score_n_hess(pf_fit, use_O_n_sq = TRUE),
                "Using '.lung' as the 'data' argument", fixed = TRUE)
  comp_obj$set_n_particles(N_fw = 250L, N_first = 250L)
  o <- comp_obj$get_get_score_n_hess()
  expect_known_value(o, "PF_get_score_n_hess-help-res-O_N_sq-old.RDS")

  # check that we get the same when we only request the score
  o_score <- comp_obj$get_get_score_n_hess(only_score = TRUE)
  expect_equal(o$score, o_score$score)
})
