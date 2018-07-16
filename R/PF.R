PF_effective_sample_size <- function(object){
  sapply(object[
    c("forward_clouds", "backward_clouds", "smoothed_clouds")],
    function(x){
    sapply(lapply(x,
      "[[", "weights"), function(z) 1 / sum(z^2))
  })
}

#' @title EM Estimation with Particle Filters and Smoothers
#' @description Method to estimate the hyper parameters with an EM algorithm.
#'
#' @inheritParams ddhazard
#' @param control see \code{\link{PF_control}}.
#' @param trace argument to get progress information. Zero will yield no info and larger integer values will yield incrementally more information.
#' @param model either \code{'logit'} for binary outcomes or \code{'exponential'} for piecewise constant exponential distributed arrival times.
#' @param seed seed to set at the start of every EM iteration.
#' @param type type of state model. Either \code{"RW"} for a [R]andom [W]alk or
#' "VAR" for [V]ector [A]uto[R]egression.
#' @param Fmat starting value for \eqn{F} when \code{type = "VAR"}. See
#' 'Details'.
#' @param ... optional way to pass arguments to \code{control}.
#'
#' @details
#' Estimates a state model of the form
#'
#' \deqn{\alpha_t = F \alpha_t + R\epsilon_t, \qquad \epsilon_t \sim N(0, Q)}
#'
#' where \eqn{F\in{\rm I\!R}^{p\times p}} has full rank,
#' \eqn{\alpha_t\in {\rm I\!R}^p}, \eqn{\epsilon_t\in {\rm I\!R}^r},
#' \eqn{r \leq p}, and \eqn{R = (e_{l_1},e_{l_2},\dots,e_{l_r})}
#' where \eqn{e_k} is column from the \eqn{p} dimensional identity matrix and
#' \eqn{l_1<l_2<\dots<l_r}. The time zero state is drawn from
#'
#' \deqn{\alpha_0\sim N(a_0, Q_0)}
#'
#' with \eqn{Q_0 \in {\rm I\!R}^{p\times p}}. The latent states,
#' \eqn{\alpha_t}, are related to the output throught the linear predictors
#'
#' \deqn{\eta_{it} = X_t(R^\top\alpha_t) + Z_t\beta}
#'
#' where \eqn{X_t\in{\rm I\!R}^{n_t\times r}} and
#' \eqn{Z_t{\rm I\!R}^{n_t\times c}} are design matrices and the outcome for a
#' individual \eqn{i} at time \eqn{t} is distributed according
#' to an exponential family member given \eqn{\eta_{it}}. \eqn{\beta} are
#' constant coefficients.
#'
#' See \code{vignette("Particle_filtering", "dynamichazard")} for details.
#'
#' @return
#' An object of class \code{PF_EM}.
#'
#' @seealso \code{\link{PF_forward_filter}} to get a more precise estimate of
#' the final log-likelihood.
#'
#' @examples
#'#####
#'# Fit model with lung data set from survival
#'# Warning: this has a longer computation time
#'
#'\dontrun{
#'library(dynamichazard)
#'.lung <- lung[!is.na(lung$ph.ecog), ]
#'set.seed(43588155)
#'pf_fit <- PF_EM(
#'  Surv(time, status == 2) ~ ph.ecog + age,
#'  data = .lung, by = 50, id = 1:nrow(.lung),
#'  Q_0 = diag(1, 3), Q = diag(1, 3),
#'  max_T = 800,
#'  control = PF_control(
#'    N_fw_n_bw = 500,
#'    N_first = 2500,
#'    N_smooth = 2500,
#'    n_max = 50,
#'    n_threads = max(parallel::detectCores(logical = FALSE), 1)),
#'  trace = 1)
#'
#'# Plot state vector estimates
#'plot(pf_fit, cov_index = 1)
#'plot(pf_fit, cov_index = 2)
#'plot(pf_fit, cov_index = 3)
#'
#'# Plot log-likelihood
#'plot(pf_fit$log_likes)
#'}
#'
#'#####
#'# Can be compared with this example from ?coxph in R 3.4.1. Though, the above
#'# only has a linear effect for age
#'
#'\dontrun{
#'cox <- coxph(
#'  Surv(time, status) ~ ph.ecog + tt(age), data= .lung,
#'  tt=function(x,t,...) pspline(x + t/365.25))
#'cox}
#'
#'
#'\dontrun{
#'######
#'# example with fixed effects
#'
#'# prepare data
#'temp <- subset(pbc, id <= 312, select=c(id, sex, time, status, edema, age))
#'pbc2 <- tmerge(temp, temp, id=id, death = event(time, status))
#'pbc2 <- tmerge(pbc2, pbcseq, id=id, albumin = tdc(day, albumin),
#'               protime = tdc(day, protime), bili = tdc(day, bili))
#'pbc2 <- pbc2[, c("id", "tstart", "tstop", "death", "sex", "edema",
#'                 "age", "albumin", "protime", "bili")]
#'pbc2 <- within(pbc2, {
#'  log_albumin <- log(albumin)
#'  log_protime <- log(protime)
#'  log_bili <- log(bili)
#'})
#'
#'# standardize
#'for(c. in c("age", "log_albumin", "log_protime", "log_bili"))
#'  pbc2[[c.]] <- drop(scale(pbc2[[c.]]))
#'
#'# fit model with extended Kalman filter
#'ddfit <- ddhazard(
#'  Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() + ddFixed(age) +
#'    ddFixed(edema) + ddFixed(log_albumin) + ddFixed(log_protime) + log_bili,
#'  pbc2, Q_0 = 100, Q = 1e-2, by = 100, id = pbc2$id,
#'  model = "exponential", max_T = 3600,
#'  control = list(eps = 1e-5, NR_eps = 1e-4, n_max = 1e4))
#'summary(ddfit)
#'
#'# fit model with particle filter
#'set.seed(88235076)
#'ppfit <- PF_EM(
#'  Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() + ddFixed(age) +
#'    ddFixed(edema) + ddFixed(log_albumin) + ddFixed(log_protime) + log_bili,
#'  pbc2, Q_0 = 100, Q = ddfit$Q * 100, # use estimate from before
#'  by = 100, id = pbc2$id,
#'  model = "exponential", max_T = 3600,
#'  control = PF_control(
#'    N_fw_n_bw = 250, N_smooth = 500, N_first = 1000, eps = 1e-3,
#'    method = "AUX_normal_approx_w_cloud_mean",
#'    n_max = 25, # just take a few iterations as an example
#'    n_threads = max(parallel::detectCores(logical = FALSE), 1)), trace = TRUE)
#'
#'# compare results
#'plot(ddfit)
#'plot(ppfit)
#'sqrt(ddfit$Q * 100)
#'sqrt(ppfit$Q)
#'rbind(ddfit$fixed_effects, ppfit$fixed_effects)
#'}
#' @export
PF_EM <- function(
  formula, data, model = "logit", by, max_T, id, a_0, Q_0, Q, order = 1,
  control = PF_control(...), trace = 0, seed = NULL, type = "RW", Fmat,
  ...){
  #####
  # checks
  if(length(order) == 1 && order != 1)
    stop(sQuote('order'), " not equal to 1 is not supported")

  if(is.character(model) && length(model) == 1 &&
     !model %in% c("logit", "exponential"))
    stop(sQuote('model'), " is not supported")

  if(missing(id)){
    if(trace > 0)
      warning("You did not parse and Id argument")
    id = 1:nrow(data)
  }

  if(is.character(type) && length(type) == 1 &&
     !type %in% c("RW", "VAR"))
    stop("Invalid ", sQuote("type"), " argument")

  if(!missing(Fmat) && type != "VAR")
    stop(sQuote("Fmat"), " should not be passed for type ", sQuote(type))

  #####
  # check if `control` has all the needed elements or if is called as in
  # version 0.5.1 or earlier
  if(!all(c(
    "N_fw_n_bw", "N_smooth", "N_first", "eps",
    "forward_backward_ESS_threshold", "method", "n_max", "n_threads",
    "smoother") %in% names(control))){
    # TODO: remove this warning some time post release 0.5.1
    warning("Please, use the ", sQuote("PF_control"), " function for the ",
            sQuote("control"), " argument.")
    control <- do.call(PF_control, control)

  }

  #####
  # find risk set and design matrix
  static_args <- .get_PF_static_args(
    formula = formula, data = data, by = by,
    max_T = if(missing(max_T)) NULL else max_T, id = id,
    trace = trace, model, order = order)

  #####
  # find matrices for state equation
  start_coefs <- get_start_values(
    formula = formula, data = data, max_T = max_T, X = static_args$X,
    fixed_terms = static_args$fixed_terms, risk_set = static_args$risk_obj,
    verbose = trace > 0, n_threads = control$n_threads, model = model,
    a_0 = if(missing(a_0)) NULL else a_0, order = order)
  a_0 <- start_coefs$a_0
  fixed_parems <- start_coefs$fixed_parems_start

  model_args <- get_state_eq_matrices(
    order = order, n_params = nrow(static_args$X),
    n_fixed = static_args$n_fixed, Q_0 = if(missing(Q_0)) NULL else Q_0,
    Q = if(missing(Q)) NULL else Q, a_0 = a_0,
    F. = if(missing(Fmat)) NULL else Fmat, type = type)
  model_args[c("L", "indicies_fix")] <- NULL

  .check_filter_input(
    Q = model_args$Q, Q_0 = model_args$Q_0, F. = model_args$F.,
    R = model_args$R, a_0 = model_args$a_0, Q_tilde = control$Q_tilde,
    fixed_parems = start_coefs$fixed_parems_start, est_fixed_in_E = FALSE,
    X = static_args$X, fixed_terms = static_args$fixed_terms, order = order)

  #####
  # build up call with symbols to get neater call stack incase of an error
  out <- .PF_EM(
    trace = trace, seed = seed, fixed_parems = fixed_parems,
    type = type, n_fixed_terms_in_state_vec =
      static_args$n_fixed_terms_in_state_vec, X = static_args$X,
    fixed_terms = static_args$fixed_terms, tstart = static_args$tstart,
    tstop = static_args$tstop, risk_obj = static_args$risk_obj,
    debug = static_args$debug, model = static_args$model, Q = model_args$Q,
    Q_0 = model_args$Q_0, F. = model_args$F., R = model_args$R,
    a_0 = model_args$a_0, N_fw_n_bw = control$N_fw_n_bw,
    N_smooth = control$N_smooth, N_first = control$N_first, eps = control$eps,
    forward_backward_ESS_threshold = control$forward_backward_ESS_threshold,
    method = control$method, n_max = control$n_max,
    n_threads = control$n_threads, smoother = control$smoother,
    Q_tilde = control$Q_tilde, est_a_0 = control$est_a_0)

  out$call <- match.call()
  out
}

#' Forward Particle Filter
#'
#' @description
#' Functions to only use the forward particle filter. Useful for log-likelihood
#' evaluation though there is an \eqn{O(d^2)} variance of the estimate where \eqn{d} is the number of time
#' periods. The number of particles specified in the \code{control} argument
#' has no effect.
#'
#' The function does not alter the \code{\link{.Random.seed}} to make sure the
#' same \code{rng.kind} is kept after the call. See \code{\link{PF_EM}} for
#' model details.
#'
#' @inheritParams PF_EM
#' @param x an \code{PF_EM} or \code{formula} object.
#' @param fixed_effects values for the fixed parameters.
#' @param N_fw number of particles.
#' @param N_first number of time zero particles to draw.
#' @param R \eqn{R} matrix in the model. See \code{\link{PF_EM}}.
#' @param Fmat \eqn{F} matrix in the model. See \code{\link{PF_EM}}.
#' @param seed \code{.GlobalEnv$.Random.seed} to set. Not \code{seed} as in
#' \code{\link{set.seed}} function. Can be used with the
#' \code{\link{.Random.seed}} returned by \code{\link{PF_EM}}.
#'
#' @return
#' An object of class \code{PF_clouds}.
#'
#' @examples
#' \dontrun{
#' # head-and-neck cancer study data. See Efron, B. (1988) doi:10.2307/2288857
#' is_censored <- c(
#'   6, 27, 34, 36, 42, 46, 48:51, 51 + c(15, 30:28, 33, 35:37, 39, 40, 42:45))
#' head_neck_cancer <- data.frame(
#'   id = 1:96,
#'   stop = c(
#'     1, 2, 2, rep(3, 6), 4, 4, rep(5, 8),
#'     rep(6, 7), 7, 8, 8, 8, 9, 9, 10, 10, 10, 11, 14, 14, 14, 15, 18, 18, 20,
#'     20, 37, 37, 38, 41, 45, 47, 47,
#'     2, 2, 3, rep(4, 4), rep(5, 5), rep(6, 5),
#'     7, 7, 7, 9, 10, 11, 12, 15, 16, 18, 18, 18, 21,
#'     21, 24, 25, 27, 36, 41, 44, 52, 54, 59, 59, 63, 67, 71, 76),
#'   event = !(1:96 %in% is_censored),
#'   group = factor(c(rep(1, 45 + 6), rep(2, 45))))
#'
#' # fit model
#' set.seed(61364778)
#' ctrl <- PF_control(
#'   N_fw_n_bw = 500, N_smooth = 2500, N_first = 2000,
#'   n_max = 1, # set to one as an example
#'   n_threads = max(parallel::detectCores(logical = FALSE), 1),
#'   eps = .001, Q_tilde = as.matrix(.3^2), est_a_0 = FALSE)
#' pf_fit <- suppressWarnings(
#'   PF_EM(
#'     survival::Surv(stop, event) ~ ddFixed(group),
#'     data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2, control = ctrl,
#'     max_T = 30))
#'
#' # the log-likelihood in the final iteration
#' (end_log_like <- tail(pf_fit$log_likes, 1))
#'
#' # gives the same
#' fw_ps <- PF_forward_filter(
#'   survival::Surv(stop, event) ~ ddFixed(group), N_fw = 500, N_first = 2000,
#'   data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2, Fmat = 1, R = 1,
#'   a_0 = pf_fit$a_0, fixed_effects = -0.5370051,
#'   control = ctrl, max_T = 30, seed = pf_fit$seed)
#' all.equal(end_log_like, logLik(fw_ps))
#'
#' # will differ since we use different number of particles
#' fw_ps <- PF_forward_filter(
#'   survival::Surv(stop, event) ~ ddFixed(group), N_fw = 1000, N_first = 3000,
#'   data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2, Fmat = 1, R = 1,
#'   a_0 = pf_fit$a_0, fixed_effects = -0.5370051,
#'   control = ctrl, max_T = 30, seed = pf_fit$seed)
#' all.equal(end_log_like, logLik(fw_ps))
#'
#' # will differ since we use the final estimates
#' fw_ps <- PF_forward_filter(pf_fit, N_fw = 500, N_first = 2000)
#' all.equal(end_log_like, logLik(fw_ps))
#' }
#' @export
PF_forward_filter <- function (x, N_fw, N_first, ...)
  UseMethod("PF_forward_filter", x)

#' @describeIn PF_forward_filter Forward particle filter with
#' \code{\link{PF_EM}} results.
#' @export
PF_forward_filter.PF_EM <- function(x, N_fw, N_first, ...){
  cl <- x$call
  cl <- cl[c(1, match(
    c("formula", "data", "model", "by", "max_T", "id", "control",
      formals(PF_control), "type", "trace", "Q_0"),
    names(cl), 0))]
  names(cl)[2] <- "x"
  xSym <- substitute(x)
  cl[ c("seed", "Fmat", "a_0", "Q", "R", "fixed_effects")] <-
    lapply(
      c("seed",    "F", "a_0", "Q", "R", "fixed_effects"),
      function(z) substitute(y$z, list(y = xSym, z = as.symbol(z))))
  cl[c("N_fw", "N_first")] <- list(N_fw, N_first)
  cl[[1]] <- quote(PF_forward_filter)

  eval(cl, parent.frame())
}

#' @describeIn PF_forward_filter Forward particle filter with formula input.
#' @export
PF_forward_filter.formula <- function(
  x, N_fw, N_first, data, model = "logit", by, max_T, id, a_0, Q_0, Q, R,
  fixed_effects, control = PF_control(...), seed = NULL, trace = 0,
  type = "RW", Fmat, ...){
  stopifnot(length(N_fw) == 1, length(N_first) == 1)

  order <- 1
  if(missing(id))
    id = 1:nrow(data)

  static_args <- .get_PF_static_args(
    formula = x, data = data, by = by,
    max_T = if(missing(max_T)) NULL else max_T, id = id,
    trace = trace, model, order = order)

  # make sure intputs are matrices if scalars are passed
  . <- function(x)
    eval(
      substitute(
        if(!is.matrix(X) && length(X) == 1) X <- as.matrix(X),
        list(X = substitute(x))),
      envir = parent.frame())
  .(Q_0)
  .(Q)
  .(R)
  .(Fmat)

  .check_filter_input(
    Q = Q, Q_0 = Q_0, F. = Fmat, R = R, a_0 = a_0, Q_tilde = control$Q_tilde,
    fixed_parems = fixed_effects, est_fixed_in_E = FALSE,
    X = static_args$X, fixed_terms = static_args$fixed_terms, order = order)
  Q_tilde <- if(is.null(control$Q_tilde))
    diag(0., ncol(Q)) else control$Q_tilde

  # set the seed
  old_seed <- .GlobalEnv$.Random.seed
  # to make sure the user has the same `rng.kind`
  on.exit(.GlobalEnv$.Random.seed <- old_seed)
  if(!is.null(seed)){
    stopifnot(length(seed) > 1) # make sure user did not use seed as in
                                # `set.seed`
    .GlobalEnv$.Random.seed <- seed
  }

  out <- particle_filter(
    fixed_parems = fixed_effects, type = type, n_fixed_terms_in_state_vec =
      static_args$n_fixed_terms_in_state_vec, X = static_args$X,
    fixed_terms = static_args$fixed_terms, tstart = static_args$tstart,
    tstop = static_args$tstop, risk_obj = static_args$risk_obj,
    debug = static_args$debug, model = static_args$model, Q = Q, Q_0 = Q_0,
    F = Fmat, R = R, is_forward = TRUE, a_0 = a_0, N_fw_n_bw = N_fw,
    N_first = N_first,
    forward_backward_ESS_threshold = control$forward_backward_ESS_threshold,
    method = control$method, n_threads = control$n_threads, Q_tilde = Q_tilde)

  structure(list(
    forward_clouds = out, backward_clouds = list(), smoothed_clouds = list(),
    transition_likelihoods = list()), class = "PF_clouds")
}

#' @importFrom graphics plot
.PF_EM <- function(
  n_fixed_terms_in_state_vec, X, fixed_terms, tstart, tstop, Q_0, Q, a_0, F.,
  R, risk_obj, n_max, n_threads, N_fw_n_bw, N_smooth, N_first, eps,
  forward_backward_ESS_threshold = NULL, debug = 0, trace,
  method = "AUX_normal_approx_w_particles", seed = NULL, smoother, model,
  fixed_parems, type, Q_tilde, est_a_0){
  cl <- match.call()
  n_vars <- nrow(X)
  fit_call <- cl
  fit_call[[1]] <- as.name("PF_smooth")

  if(is.null(Q_tilde))
    fit_call[["Q_tilde"]] <- diag(0, n_vars)
  fit_call[["F"]]   <- eval(fit_call[["F."]] , parent.frame())
  fit_call[["a_0"]] <- eval(fit_call[["a_0"]], parent.frame())
  fit_call[["Q"]]   <- eval(fit_call[["Q"]]  , parent.frame())

  fit_call[c("eps", "seed", "F.", "trace", "est_a_0")] <- NULL

  #####
  # set the seed as in r-source/src/library/stats/R/lm.R `simulate.lm`
  if(is.null(seed)){
    if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      runif(1)                     # initialize the RNG if necessary
    seed <- get(".Random.seed", envir = .GlobalEnv)

  } else {
    set.seed(seed)
    seed <- get(".Random.seed", envir = .GlobalEnv)

  }

  log_likes <- rep(NA_real_, n_max)
  log_like <- log_like_max <- -Inf
  for(i in 1:n_max){
    if(trace > 0){
      if(i != 1)
        cat("\n\n\n")
      cat("#######################\nStarting EM iteration", i, "\n")

      cat("a_0 is:\n")
      print(fit_call$a_0)

      if(length(fixed_parems) > 0){
        cat("Fixed parameters are:\n")
        print(fixed_parems)
      }
      if(type == "VAR"){
        cat("F is:\n")
        print(fit_call$F)
      }

      cat("chol(Q) is:\n")
      print(chol(fit_call$Q))
    }
    log_like_old <- log_like
    log_like_max <- max(log_like, log_like_max)

    #####
    # find clouds
    assign(".Random.seed", seed, envir = .GlobalEnv)
    clouds <- eval(fit_call, envir = parent.frame())

    if(trace > 0){
      cat("Plotting state vector mean and quantiles for iteration", i, "\n")
      plot(clouds, main = paste0("EM iteration ", i))
      plot(clouds, type = "forward_clouds", add = TRUE, qlvls = c(), lty = 2)
      plot(clouds, type = "backward_clouds", add = TRUE, qlvls = c(), lty = 3)

      cat("Effective sample sizes are:\n")
      print(effective_sample_size <- PF_effective_sample_size(clouds))
    }

    #####
    # update parameters in state equation
    if(trace > 0)
      cat("Updating parameters in state model...\n")
    a_0_old <- fit_call$a_0
    Q_old <- fit_call$Q
    F_old <- fit_call$F

    if(type == "RW"){
      sum_stats <- compute_summary_stats_first_o_RW(
        clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
        debug = trace > 2)
      if(est_a_0)
        a_0 <- drop(sum_stats[[1]]$E_xs)
      Q <- Reduce("+", lapply(sum_stats, "[[", "E_x_less_x_less_one_outers"))
      Q <- Q / length(sum_stats)

      fit_call$a_0 <- a_0
      fit_call$Q <- Q
    } else if (type == "VAR") {
      new_params <- PF_est_params_dens(
        clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
        debug = trace > 1)
      fit_call$F <- F. <- new_params$R_top_F # TODO: need to change for higher
                                             #       order models
      fit_call$Q <- Q <- new_params$Q

    } else
      stop(sQuote("type"), " not implemented")

    #####
    # Update fixed effects
    if(length(fixed_parems) > 0){
      if(trace > 0)
        cat("Updating fixed effects...\n")

      fit_call$fixed_parems <- fixed_parems <- .PF_update_fixed(
        clouds = clouds$smoothed_clouds, risk_obj = risk_obj, model = model,
        R = R, X = X, fixed_terms = fixed_terms, fixed_parems = fixed_parems,
        nthreads = n_threads, tstart = tstart, tstop = tstop,
        debug = trace > 1L)
    }

    #####
    # compute log likelihood and check for convergernce
    log_likes[i] <- log_like <- logLik(clouds)

    if(trace > 0)
      cat("The log likelihood in iteration ", i, " is ", log_like,
          ". Largest log likelihood before this iteration is ", log_like_max,
          "\n", sep = "")

    Q_relative_norm <- norm(Q_old - Q) / (norm(Q_old) + 1e-8)
    a_0_relative_norm <- norm(t(a_0 - a_0_old)) / (norm(t(a_0_old)) + 1e-8)
    F_norm <- norm(F_old - F.) / (norm(F_old) + 1e-8)

    if(trace > 0){
      msg <- "The relative norm of the change in"
      if(type == "VAR")
        msg <- paste(msg, "F,")

      msg <- paste(msg, "a_0 and Q are")
      if(type == "VAR")
        msg <- paste0(msg, " ", sprintf("%.5f", F_norm), ",")

      msg <- paste(
        msg, sprintf("%.5f", a_0_relative_norm), "and",
        sprintf("%.5f", Q_relative_norm), "at iteration", i)
      cat(msg)
    }

    if(has_converged <-
       Q_relative_norm < eps &&
       a_0_relative_norm < eps &&
       (type != "VAR" || F_norm < eps))
      break
  }

  if(!has_converged)
    warning("Method did not converge.")

  if(!exists("effective_sample_size", envir = environment()))
    effective_sample_size <- PF_effective_sample_size(clouds)

  return(structure(list(
    call = cl, clouds = clouds, a_0 = a_0, fixed_effects = fixed_parems, Q = Q,
    F = fit_call$F, R = R,
    summary_stats = if(type == "RW") sum_stats else NULL,
    log_likes = log_likes[1:i], n_iter = i,
    effective_sample_size = effective_sample_size, seed = seed),
    class = "PF_EM"))
}

#' @title Auxiliary for Controlling Particle Fitting
#'
#' @description
#' Auxiliary for additional settings with \code{\link{PF_EM}}.
#'
#' @param N_fw_n_bw number of particles to use in forward and backward filter.
#' @param N_smooth number of particles to use in particle smoother.
#' @param N_first number of particles to use at time \eqn{0} and time \eqn{d + 1}.
#' @param eps convergence threshold in EM method.
#' @param forward_backward_ESS_threshold required effective sample size to not re-sample in the particle filters.
#' @param method method for forward, backward and smoothing filter.
#' @param n_max maximum number of iterations of the EM algorithm.
#' @param n_threads maximum number threads to use in the computations.
#' @param smoother smoother to use.
#' @param Q_tilde covariance matrix of additional error term to add to the
#' proposal distributions. \code{NULL} implies no additional error term.
#' @param est_a_0 \code{FALSE} if the starting value of the state model should
#' be fixed. Does not apply for \code{type = "VAR"}.
#'
#' @return
#' A list with components named as the arguments.
#'
#' @seealso
#' \code{\link{PF_EM}}
#'
#' @export
PF_control <- function(
  N_fw_n_bw = NULL, N_smooth = NULL, N_first = NULL,
  eps = 1e-2, forward_backward_ESS_threshold = NULL,
  method = "AUX_normal_approx_w_cloud_mean", n_max = 25,
  n_threads = getOption("ddhazard_max_threads"), smoother = "Fearnhead_O_N",
  Q_tilde = NULL, est_a_0 = TRUE){
  control <- list(
    N_fw_n_bw = N_fw_n_bw, N_smooth = N_smooth, N_first = N_first, eps = eps,
    forward_backward_ESS_threshold = forward_backward_ESS_threshold,
    method = method, n_max = n_max, n_threads = n_threads, smoother = smoother,
    Q_tilde = Q_tilde, est_a_0 = est_a_0)

  check_n_particles_expr <- function(N_xyz)
    eval(bquote({
      if(is.null(control[[.(N_xyz)]]))
        stop("Please supply the number of particle in ", sQuote(paste0(
          "PF_control(", .(N_xyz), ")")))
    }), parent.frame())

  check_n_particles_expr("N_first")
  check_n_particles_expr("N_fw_n_bw")
  check_n_particles_expr("N_smooth")

  return(control)
}


.get_PF_static_args <- function(
  formula, data, by, max_T = NULL, id, trace, model, order){
  # get design matrix and risk set
  tmp <- get_design_matrix_and_risk_obj(
    formula = formula, data = data, by = by,
    max_T = if(is.null(max_T)) NULL else max_T, verbose = trace > 0,
    is_for_discrete_model = model == "logit", id = id)

  if(trace > 0)
    report_pre_liminary_stats_before_EM(
      risk_set = tmp$risk_set, Y = tmp$X_Y$Y)

  # tranpose due to column-major storage and we want to look up individuals
  tmp$X_Y$X           <- t(tmp$X_Y$X)
  tmp$X_Y$fixed_terms <- t(tmp$X_Y$fixed_terms)

  with(tmp, list(
    n_fixed_terms_in_state_vec = 0, X = X_Y$X, fixed_terms = X_Y$fixed_terms,
    tstart = X_Y$Y[, 1], tstop = X_Y$Y[, 2], risk_obj = risk_set,
    debug = max(0, trace - 1), model = model))
}


.PF_update_fixed <- function(
  clouds, risk_obj, R, X, fixed_terms, fixed_parems, model, nthreads,
  tstart, tstop, debug){
  if(!model %in% c("logit", "exponential"))
    stop(sQuote(model), " is not implemented with fixed effects")

  family_arg <- switch(
    model,
    logit = "binomial",
    exponential = "poisson")

  out <- NULL
  R_top <- t(R)
  for(i in 1:length(clouds)){
    cl <- clouds[[i]]
    risk_set <- risk_obj$risk_sets[[i]]

    ran_vars <- X[, risk_set, drop = FALSE]
    X_i <- fixed_terms[, risk_set, drop = FALSE]
    y_i <- risk_obj$is_event_in[risk_set] == (i - 1)

    good <- which(drop(cl$weights) >= 1e-7)
    dts <-
      pmin(tstop[risk_set], risk_obj$event_times[i + 1]) -
      pmax(tstart[risk_set], risk_obj$event_times[i])

    ws <- cl$weights[good] # TODO: maybe scale up the weights at this point?
    particle_coefs <- R_top %*% cl$states[, good, drop = FALSE]

    out <- c(out, list(
      pf_fixed_effect_iteration(
        X = X_i, Y = y_i, dts = dts, cloud = particle_coefs,
        cl_weights = ws, ran_vars = ran_vars, beta = fixed_parems,
        family = family_arg, max_threads = nthreads, debug = debug)))
  }

  f_stack <- do.call(c, lapply(out, "[[", "f"))
  R_stack <- lapply(out, .get_R)
  R_stack <- do.call(rbind, R_stack)

  qr. <- qr(R_stack, LAPACK = TRUE)
  f <- qr.qty(qr., f_stack)[1:nrow(X_i)]

  out <- list(list(
    R = qr.R(qr.), f = f,
    pivot = qr.$pivot - 1)) # less one to have zero index as cpp code

  R <- .get_R(out[[1]])
  drop(solve(t(R) %*% R, t(R) %*% out[[1]]$f))
}

.get_R <- function(o){
  piv <- drop(o$pivot) + 1
  piv[piv] <- 1:length(piv)
  o$R[, piv, drop = FALSE]
}
