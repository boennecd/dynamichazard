utils::globalVariables("is_restricted")

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
#' @param model either \code{'logit'} for binary outcomes with the logistic
#' link function, \code{'cloglog'} for binary outcomes with the inverse cloglog
#' link function, or \code{'exponential'} for piecewise constant exponential distributed arrival times.
#' @param seed seed to set at the start of every EM iteration. See
#' \code{\link{set.seed}}.
#' @param type type of state model. Either \code{"RW"} for a [R]andom [W]alk or
#' "VAR" for [V]ector [A]uto[R]egression.
#' @param Fmat starting value for \eqn{F} when \code{type = "VAR"}. See
#' 'Details' in \code{\link{PF_EM}}.
#' @param fixed_effects starting values for fixed effects if any. See
#' \code{\link{ddFixed}}.
#' @param G,theta,J,K,psi,phi parameters for a restricted \code{type = "VAR"} model.
#' See the vignette mentioned in 'Details' of \code{\link{PF_EM}} and the
#' examples linked to in 'See Also'.
#' @param fixed  two-sided \code{\link{formula}} to be used
#' with \code{random} instead of \code{formula}. It is of the form
#' \code{Surv(tstart, tstop, event) ~ x} or
#' \code{Surv(tstart, tstop, event) ~ - 1} for no fixed effects.
#' @param random one-sided \code{\link{formula}} to be used
#' with \code{fixed} instead of \code{formula}. It is of the form
#' \code{~ z}.
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
#' \eqn{\alpha_t}, are related to the output through the linear predictors
#'
#' \deqn{\eta_{it} = X_t(R^+\alpha_t) + Z_t\beta}
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
#' See the examples at https://github.com/boennecd/dynamichazard/tree/master/examples.
#'
#' @section Warning:
#' The function is still under development so the output and API may change.
#'
#' @examples
#'\dontrun{
#'#####
#'# Fit model with lung data set from survival
#'# Warning: long-ish computation time
#'
#'library(dynamichazard)
#' .lung <- lung[!is.na(lung$ph.ecog), ]
#' # standardize
#' .lung$age <- scale(.lung$age)
#'
#' # fit
#' set.seed(43588155)
#' pf_fit <- PF_EM(
#'  Surv(time, status == 2) ~ ddFixed(ph.ecog) + age,
#'  data = .lung, by = 50, id = 1:nrow(.lung),
#'  Q_0 = diag(1, 2), Q = diag(.5^2, 2),
#'  max_T = 800,
#'  control = PF_control(
#'     N_fw_n_bw = 500, N_first = 2500, N_smooth = 5000,
#'     n_max = 50, eps = .001, Q_tilde = diag(.2^2, 2), est_a_0 = FALSE,
#'     n_threads = max(parallel::detectCores(logical = FALSE), 1)))
#'
#'# Plot state vector estimates
#'plot(pf_fit, cov_index = 1)
#'plot(pf_fit, cov_index = 2)
#'
#'# Plot log-likelihood
#'plot(pf_fit$log_likes)
#'}
#'\dontrun{
#'######
#'# example with fixed intercept
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
#'  control = ddhazard_control(eps = 1e-5, NR_eps = 1e-4, n_max = 1e4))
#'summary(ddfit)
#'
#'# fit model with particle filter
#'set.seed(88235076)
#'pf_fit <- PF_EM(
#'   Surv(tstart, tstop, death == 2) ~ ddFixed_intercept() + ddFixed(age) +
#'     ddFixed(edema) + ddFixed(log_albumin) + ddFixed(log_protime) + log_bili,
#'   pbc2, Q_0 = 2^2, Q = ddfit$Q * 100, # use estimate from before
#'   by = 100, id = pbc2$id,
#'   model = "exponential", max_T = 3600,
#'   control = PF_control(
#'     N_fw_n_bw = 500, N_smooth = 2500, N_first = 1000, eps = 1e-3,
#'     method = "AUX_normal_approx_w_cloud_mean", est_a_0 = FALSE,
#'     Q_tilde = as.matrix(.1^2),
#'     n_max = 25, # just take a few iterations as an example
#'     n_threads = max(parallel::detectCores(logical = FALSE), 1)))
#'
#'# compare results
#'plot(ddfit)
#'plot(pf_fit)
#'sqrt(ddfit$Q * 100)
#'sqrt(pf_fit$Q)
#'rbind(ddfit$fixed_effects, pf_fit$fixed_effects)
#'}
#'\dontrun{
#' #####
#' # simulation example with `random` and `fixed` argument and a restricted
#' # model
#'
#' # g groups with k individuals in each
#' g <- 3L
#' k <- 400L
#'
#' # matrices for state equation
#' p <- g + 1L
#' G <- matrix(0., p^2, 2L)
#' for(i in 1:p)
#'   G[i + (i - 1L) * p, 1L + (i == p)] <- 1L
#'
#' theta <- c(.9, .8)
#' # coefficients in transition density
#' (F. <- matrix(as.vector(G %*% theta), 4L, 4L))
#'
#' J <- matrix(0., ncol = 2L, nrow = p)
#' J[-p, 1L] <- J[p, 2L] <- 1
#' psi <- c(log(c(.3, .1)))
#'
#' K <- matrix(0., p * (p - 1L) / 2L, 2L)
#' j <- 0L
#' for(i in (p - 1L):1L){
#'   j <- j + i
#'   K[j, 2L] <- 1
#' }
#' K[K[, 2L] < 1, 1L] <- 1
#' phi <- log(-(c(.8, .3) + 1) / (c(.8, .3) - 1))
#'
#' V <- diag(exp(drop(J %*% psi)))
#' C <- diag(1, ncol(V))
#' C[lower.tri(C)] <- 2/(1 + exp(-drop(K %*% phi))) - 1
#' C[upper.tri(C)] <- t(C)[upper.tri(C)]
#' (Q <- V %*% C %*% V)     # covariance matrix in transition density
#' cov2cor(Q)
#'
#' Q_0 <- get_Q_0(Q, F.)    # time-invariant covariance matrix
#' beta <- c(rep(-6, g), 0) # all groups have the same long run mean intercept
#'
#' # simulate state variables
#' set.seed(56219373)
#' n_periods <- 300L
#' alphas <- matrix(nrow = n_periods + 1L, ncol = p)
#' alphas[1L, ] <- rnorm(p) %*% chol(Q_0)
#' for(i in 1:n_periods + 1L)
#'   alphas[i, ] <- F. %*% alphas[i - 1L, ] + drop(rnorm(p) %*% chol(Q))
#'
#' alphas <- t(t(alphas) + beta)
#'
#' # plot state variables
#' matplot(alphas, type = "l", lty = 1)
#'
#' # simulate individuals' outcome
#' n_obs <- g * k
#' df <- lapply(1:n_obs, function(i){
#'   # find the group
#'   grp <- (i - 1L) %/% (n_obs / g) + 1L
#'
#'   # left-censoring
#'   tstart <- max(0L, sample.int((n_periods - 1L) * 2L, 1) - n_periods + 1L)
#'
#'   # covariates
#'   x <- c(1, rnorm(1))
#'
#'   # outcome (stop time and event indicator)
#'   osa <- NULL
#'   oso <- NULL
#'   osx <- NULL
#'   y <- FALSE
#'   for(tstop in (tstart + 1L):n_periods){
#'     sigmoid <- 1 / (1 + exp(- drop(x %*% alphas[tstop + 1L, c(grp, p)])))
#'     if(sigmoid > runif(1)){
#'       y <- TRUE
#'       break
#'     }
#'     if(.01 > runif(1L) && tstop < n_periods){
#'       # sample new covariate
#'       osa <- c(osa, tstart)
#'       tstart <- tstop
#'       oso <- c(oso, tstop)
#'       osx <- c(osx, x[2])
#'       x[2] <- rnorm(1)
#'     }
#'   }
#'
#'   cbind(
#'     tstart = c(osa, tstart), tstop = c(oso, tstop),
#'     x = c(osx, x[2]), y = c(rep(FALSE, length(osa)), y), grp = grp,
#'     id = i)
#' })
#' df <- data.frame(do.call(rbind, df))
#' df$grp <- factor(df$grp)
#'
#' # fit model. Start with "cheap" iterations
#' fit <- PF_EM(
#'   fixed = Surv(tstart, tstop, y) ~ x, random = ~ grp + x - 1,
#'   data = df, model = "logit", by = 1L, max_T = max(df$tstop),
#'   Q_0 = diag(1.5^2, p), id = df$id, type = "VAR",
#'   G = G, theta = c(.5, .5), J = J, psi = log(c(.1, .1)),
#'   K = K, phi = log(-(c(.4, 0) + 1) / (c(.4, 0) - 1)),
#'   control = PF_control(
#'     N_fw_n_bw = 100L, N_smooth = 100L, N_first = 500L,
#'     method = "AUX_normal_approx_w_cloud_mean",
#'     nu = 5L, # sample from multivariate t-distribution
#'     n_max = 100L,  averaging_start = 50L,
#'     smoother = "Fearnhead_O_N", eps = 1e-4, covar_fac = 1.2,
#'     n_threads = 4L # depends on your cpu(s)
#'   ),
#'   trace = 1L)
#' plot(fit$log_likes) # log-likelihood approximation at each iterations
#'
#' # take more iterations with more particles
#' cl <- fit$call
#' ctrl <- cl[["control"]]
#' ctrl[c("N_fw_n_bw", "N_smooth", "N_first", "n_max",
#'        "averaging_start")] <- list(500L, 2000L, 5000L, 200L, 30L)
#' cl[["control"]] <- ctrl
#' cl[c("phi", "psi", "theta")] <- list(fit$phi, fit$psi, fit$theta)
#' fit_extra <- eval(cl)
#'
#' plot(fit_extra$log_likes) # log-likelihood approximation at each iteration
#'
#' # check estimates
#' sqrt(diag(fit_extra$Q))
#' sqrt(diag(Q))
#' cov2cor(fit_extra$Q)
#' cov2cor(Q)
#' fit_extra$F
#' F.
#'
#' # plot predicted state variables
#' for(i in 1:p){
#'   plot(fit_extra, cov_index = i)
#'   abline(h = 0, lty = 2)
#'   lines(1:nrow(alphas) - 1, alphas[, i] - beta[i], lty = 3)
#' }
#'}
#' @export
PF_EM <- function(
  formula, data, model = "logit", by, max_T, id, a_0, Q_0, Q, order = 1,
  control = PF_control(...), trace = 0, seed = NULL, type = "RW",
  fixed = NULL, random = NULL,
  Fmat, fixed_effects, G, theta, J, K, psi, phi, ...){
  #####
  # checks
  eval(.PF_check_args_expre, environment())

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
    trace = trace, model, order = order, fixed = fixed, random = random)

  if(type == "RW" && !missing(fixed) && !missing(random))
    local({
      func <- function(x){
        v <- attr(x, "term.labels")
        if(attr(x, "intercept"))
          v <- c("intercept", v)
        v
      }

      rvars <- func(static_args$terms$random)
      fvars <- func(static_args$terms$fixed)

      un <- intersect(rvars, fvars)
      if(length(un) > 0L)
        warning(paste0(
          "The following terms are in both ", sQuote("fixed"), " and ",
          sQuote("random"), " with ", sQuote("type = \"RW\""), ":  ",
          paste0(sQuote(un), collapse = ", ")))
    })

  #####
  # find matrices for state equation
  start_coefs <- get_start_values(
    data = data, formula = formula, max_T = max_T, X = static_args$X,
    fixed_terms = static_args$fixed_terms, risk_set = static_args$risk_obj,
    verbose = trace > 0, n_threads = control$n_threads, model = model,
    a_0 = if(missing(a_0)) NULL else a_0, order = order,
    fixed_parems_start = if(missing(fixed_effects)) NULL else fixed_effects,
    fixed = fixed, random = random, type = type)
  a_0 <- start_coefs$a_0
  fixed_params <- start_coefs$fixed_parems_start

  if(is_restricted){
    model_args <- list(
      G = G, theta = theta, J = J, K = K, psi = psi, phi = phi,
      R = diag(nrow(static_args$X)), a_0 = a_0, Q_0 = Q_0)

  } else {
    model_args <- get_state_eq_matrices(
      order = order, n_params = nrow(static_args$X),
      n_fixed = static_args$n_fixed, Q_0 = if(missing(Q_0)) NULL else Q_0,
      Q = if(missing(Q)) NULL else Q, a_0 = a_0,
      F. = if(missing(Fmat)) NULL else Fmat, type = type)
    model_args[c("L", "indicies_fix")] <- NULL

  }

  .check_filter_input(
    Q = model_args$Q, Q_0 = model_args$Q_0, F. = model_args$F.,
    R = model_args$R, a_0 = model_args$a_0, Q_tilde = control$Q_tilde,
    fixed_parems = start_coefs$fixed_parems_start, est_fixed_in_E = FALSE,
    X = static_args$X, fixed_terms = static_args$fixed_terms, order = order,
    G = model_args$G, theta = model_args$theta, J = model_args$J,
    K = model_args$K, psi = model_args$psi, phi = model_args$phi)

  #####
  # build up call with symbols to get neater call stack incase of an error
  out <- .PF_EM(
    trace = trace, seed = seed, fixed_params = fixed_params,
    type = type, n_fixed_terms_in_state_vec =
      static_args$n_fixed_terms_in_state_vec, X = static_args$X,
    fixed_terms = static_args$fixed_terms, tstart = static_args$tstart,
    tstop = static_args$tstop, risk_obj = static_args$risk_obj,
    debug = static_args$debug, model = static_args$model, Q = model_args$Q,
    Q_0 = model_args$Q_0, F. = model_args$F., R = model_args$R,
    a_0 = model_args$a_0, G = model_args$G, J = model_args$J, K = model_args$K,
    theta = model_args$theta, psi = model_args$psi, phi = model_args$phi,
    N_fw_n_bw = control$N_fw_n_bw, nu = control$nu,
    N_smooth = control$N_smooth, N_smooth_final = control$N_smooth_final,
    N_first = control$N_first, eps = control$eps,
    forward_backward_ESS_threshold = control$forward_backward_ESS_threshold,
    method = control$method, n_max = control$n_max,
    n_threads = control$n_threads, smoother = control$smoother,
    Q_tilde = control$Q_tilde, est_a_0 = control$est_a_0,
    covar_fac = control$covar_fac, ftol_rel = control$ftol_rel,
    averaging_start = control$averaging_start, fix_seed = control$fix_seed)

  out <- .set_PF_names(out, rng_names = row.names(static_args$X),
                       fixed_names = rownames(static_args$fixed_terms))

  out[c("call", "terms", "has_fixed_intercept", "xlev", "control", "type")] <-
    list(
      match.call(), static_args$terms, static_args$has_fixed_intercept,
      static_args$xlev, control, type)
  if(!is.null(fixed))
    out$fixed <- fixed
  if(!is.null(random))
    out$random <- random
  out
}

.PF_check_args_expre <- expression({
  if(length(order) == 1 && order != 1)
    stop(sQuote('order'), " not equal to 1 is not supported")

  if(is.character(model) && length(model) == 1 &&
     !model %in% c("logit", "cloglog", "exponential"))
    stop(sQuote('model'), " is not supported")

  if(missing(id)){
    if(trace > 0)
      warning("You did not parse and Id argument")
    id <- 1:nrow(data)
  }

  if(is.character(type) && length(type) == 1 &&
     !type %in% c("RW", "VAR"))
    stop("Invalid ", sQuote("type"), " argument")

  if(!missing(Fmat) && type != "VAR")
    stop(sQuote("Fmat"), " should not be passed for type ", sQuote(type))

  has_restrict <- !c(missing(G), missing(theta), missing(J), missing(K),
                     missing(psi), missing(phi))
  is_restricted <- all(has_restrict)
  str_if_err <- paste0(
    sQuote("G"), ", ", sQuote("theta"), ", ", sQuote("J"), ", ", sQuote("K"),
    ", ", sQuote("psi"), ", and ", sQuote("phi"))
  if(!(missing(Q) && missing(Fmat)) && !all(!has_restrict))
    stop("Either supply ", sQuote("Q"), " and ", sQuote("Fmat"),
         " or ", str_if_err)
  if(any(has_restrict) && type != "VAR")
    stop(str_if_err, " supplied with ", sQuote("type"),
         " ", sQuote(type))
  if(any(has_restrict) && !all(has_restrict))
    stop("Missing one of ", str_if_err)

  if(is.null(fixed) != is.null(random))
    stop("supply either both ", sQuote("fixed"), " and ", sQuote("random"),
         " or none of them")

  if(!missing(formula) & !is.null(fixed) &  !is.null(random))
    stop("Use either ", sQuote("formula"), " or ", sQuote("fixed"), " and ",
         sQuote("random"))
})

.set_PF_names <- function(obj, fixed_names, rng_names){
  fn <- length(fixed_names)
  rn <- length(rng_names)

  # we assume that there is atleast one random effect
  colnames(obj$R) <- rng_names
  state_names <- character(nrow(obj$R))
  for(i in seq_len(rn)){
    j <- which(obj$R[, i] == 1)
    state_names[j] <- rng_names[i]
  }
  rownames(obj$R) <- state_names
  if(!is.null(obj$F))
    dimnames(obj$F) <- list(state_names, state_names)
  if(!is.null(obj$Q))
    dimnames(obj$Q) <- list(rng_names, rng_names)
  if(!is.null(obj$a_0))
    names(obj$a_0) <- state_names
  if(length(fixed_names) > 0)
    names(obj$fixed_effects) <- fixed_names

  if(length(obj$EM_ests) > 0){
    n_iter <- nrow(obj$EM_ests$a_0)
    em_names <- paste0("EM it ", 1:n_iter)
    dimnames(obj$EM_ests$a_0) <- list(em_names, state_names)
    if(length(fixed_names) > 0)
      dimnames(obj$EM_ests$fixed_effects) <- list(em_names, fixed_names)
    dimnames(obj$EM_ests$F) <- list(state_names, state_names, em_names)
    dimnames(obj$EM_ests$Q) <- list(rng_names, rng_names, em_names)

  }

  if(!is.null(obj$G) && !is.null(obj$theta)){
    colnames(obj$G) <- names(obj$theta) <-
      paste0("theta", seq_along(obj$theta))
    tmp <- outer(
      state_names, state_names, function(x, y) paste0(x, ":", y))
    rownames(obj$G) <- as.vector(tmp)

  }

  if(!is.null(obj$J) && !is.null(obj$psi)){
    colnames(obj$J) <- names(obj$psi) <- paste0("psi", seq_along(obj$psi))
    rownames(obj$J) <- rng_names

  }

  if(!is.null(obj$K) && !is.null(obj$phi)){
    if(length(obj$phi) > 0)
      colnames(obj$K) <- names(obj$phi) <- paste0("phi", seq_along(obj$phi))
    tmp <- outer(rng_names, rng_names, function(x, y) paste0(x, ":", y))
    rownames(obj$K) <- tmp[lower.tri(tmp)]

  }

  obj
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
#' @param seed \code{.GlobalEnv$.Random.seed} to set. Not \code{seed} as in
#' \code{\link{set.seed}} function. Can be used with the
#' \code{\link{.Random.seed}} returned by \code{\link{PF_EM}}.
#' @param G,theta,J,K,psi,phi parameters for a restricted \code{type = "VAR"} model.
#' See the vignette mentioned in 'Details' of \code{\link{PF_EM}} and the
#' examples linked to in 'See Also'.
#'
#' @return
#' An object of class \code{PF_clouds}.
#'
#' @section Warning:
#' The function is still under development so the output and API may change.
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
#' (end_log_like <- logLik(pf_fit))
#'
#' # gives the same
#' fw_ps <- PF_forward_filter(
#'   survival::Surv(stop, event) ~ ddFixed(group), N_fw = 500, N_first = 2000,
#'   data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2,
#'   a_0 = pf_fit$a_0, fixed_effects = -0.5370051,
#'   control = ctrl, max_T = 30, seed = pf_fit$seed)
#' all.equal(c(end_log_like), c(logLik(fw_ps)))
#'
#' # will differ since we use different number of particles
#' fw_ps <- PF_forward_filter(
#'   survival::Surv(stop, event) ~ ddFixed(group), N_fw = 1000, N_first = 3000,
#'   data = head_neck_cancer, by = 1, Q_0 = 1, Q = 0.1^2,
#'   a_0 = pf_fit$a_0, fixed_effects = -0.5370051,
#'   control = ctrl, max_T = 30, seed = pf_fit$seed)
#' all.equal(c(end_log_like), c(logLik(fw_ps)))
#'
#' # will differ since we use the final estimates
#' fw_ps <- PF_forward_filter(pf_fit, N_fw = 500, N_first = 2000)
#' all.equal(c(end_log_like), c(logLik(fw_ps)))
#' }
#' @export
PF_forward_filter <- function (x, N_fw, N_first, ...)
  UseMethod("PF_forward_filter", x)

#' @describeIn PF_forward_filter Forward particle filter using the
#' estimates of an \code{\link{PF_EM}} call.
#' @export
PF_forward_filter.PF_EM <- function(x, N_fw, N_first, seed, ...){
  cl <- x$call
  ma <- c(1L, match(
    c("formula", "fixed", "random", "data", "model", "by", "max_T", "id",
      "control", formals(PF_control), "type", "trace", "Q_0", "G", "J", "K"),
    names(cl), 0))
  cl <- cl[ma]
  if(ma[2L] == 0L){ # no forumla
    names(cl)[4L] <- "x"
  } else # the data argument
      names(cl)[2L] <- "x" # the formula argument
  xSym <- substitute(x)
  vas <- rbind(
    c("seed", "Fmat", "a_0", "Q", "R", "fixed_effects", "psi", "phi", "theta"),
    c("seed", "F", "a_0", "Q", "R", "fixed_effects", "psi", "phi", "theta"))
  for(i in 1:ncol(vas)){
    if(!is.null(x[[vas[2, i]]]))
      cl[[vas[1, i]]] <- substitute(
        x$y, list(x = xSym, y = as.symbol(vas[2, i])))
  }
  if(!missing(seed))
    cl[["seed"]] <- substitute(seed)
  else if(!x$control$fix_seed)
    cl[["seed"]] <- NULL

  cl[c("N_fw", "N_first")] <- list(N_fw, N_first)
  cl[[1L]] <- quote(PF_forward_filter)

  if(is.null(cl$type) || # RW is the default
     cl$type == "RW")
    cl["Fmat"] <- NULL

  if(!is.null(cl$theta))
    cl[c("Q", "Fmat")] <- NULL

  eval(cl, parent.frame())
}

#' @describeIn PF_forward_filter Forward particle filter with formula input.
#' @export
PF_forward_filter.formula <- function(
  x, N_fw, N_first, data, model = "logit", by, max_T, id, a_0, Q_0, Q,
  fixed_effects, control = PF_control(...), seed = NULL, trace = 0,
  G, theta, J, K, psi, phi, type = "RW", Fmat, ...){
  cl <- match.call()
  idx <- match(c("x", "data"), names(cl), NA_integer_)
  stopifnot(!anyNA(idx))
  names(cl)[idx] <- c("formula", "x")
  cl[[1L]] <- quote(PF_forward_filter)

  eval(cl, parent.frame())
}

#' @describeIn PF_forward_filter Forward particle filter with \code{data.frame}
#' data input as \code{x} instead of \code{data}. Can be used with \code{fixed}
#' and \code{random} argument.
#' @export
PF_forward_filter.data.frame <- function(
  x, N_fw, N_first, formula, model = "logit", by, max_T, id, a_0,
  Q_0, Q, fixed_effects, control = PF_control(...), seed = NULL, trace = 0,
  fixed = NULL, random = NULL, G, theta, J, K, psi, phi, type = "RW", Fmat,
  order = 1, ...){
  data <- x

  eval(.PF_check_args_expre, environment())
  stopifnot(length(N_fw) == 1, length(N_first) == 1)

  cl <- quote(.get_PF_static_args(
    data = data, by = by,
    max_T = if(missing(max_T)) NULL else max_T, id = id,
    trace = trace, model = model, order = order))
  if(!missing(formula))
    cl[["formula"]] <- quote(formula) else
      cl[c("fixed", "random")] <- list(quote(fixed), quote(random))

  static_args <- eval(cl, environment())

  if(type == "RW"){
    . <- function(x)
      eval(
        substitute(
          if(!is.matrix(X) && length(X) == 1) X <- as.matrix(X),
          list(X = substitute(x))),
        envir = parent.frame())
    .(Q_0)
    .(Q)
    if(order == 1L)
      Fmat <- diag(ncol(Q))

  } else if(type == "VAR"){
    if(is_restricted){
      Fmat <- .get_F(G, theta)
      Q <- .get_Q(J = J, K = K, psi = psi, phi = phi)$Q
    }

    Q_0 <- get_Q_0(Fmat = Fmat, Qmat = Q)
  }

  # set R
  if(order != 1)
    stop("Method not implemented with ", sQuote("order"), " ", order) else
      R <- diag(ncol(Q))

  .check_filter_input(
    Q = Q, Q_0 = Q_0, F. = Fmat, R = R, a_0 = a_0, Q_tilde = control$Q_tilde,
    fixed_parems = fixed_effects, est_fixed_in_E = FALSE,
    X = static_args$X, fixed_terms = static_args$fixed_terms, order = order)

  Q_tilde <- get_Q_tilde(control$Q_tilde, ncol(Q))

  if(!is.null(seed)){
    # set the seed
    old_seed <- .GlobalEnv$.Random.seed
    # to make sure the user has the same `rng.kind`
    on.exit(.GlobalEnv$.Random.seed <- old_seed)

    stopifnot(length(seed) > 1) # make sure user did not use seed as in
    # `set.seed`
    if(control$fix_seed)
      .GlobalEnv$.Random.seed <- seed
  }

  out <- particle_filter(
    fixed_params = fixed_effects, type = type, n_fixed_terms_in_state_vec =
      static_args$n_fixed_terms_in_state_vec, X = static_args$X,
    fixed_terms = static_args$fixed_terms, tstart = static_args$tstart,
    tstop = static_args$tstop, risk_obj = static_args$risk_obj,
    debug = static_args$debug, model = static_args$model, Q = Q, Q_0 = Q_0,
    F = Fmat, R = R, is_forward = TRUE, a_0 = a_0, N_fw_n_bw = N_fw,
    N_first = N_first, nu = if(is.null(control$nu)) 0L else control$nu,
    forward_backward_ESS_threshold = control$forward_backward_ESS_threshold,
    method = control$method, n_threads = control$n_threads, Q_tilde = Q_tilde,
    covar_fac = control$covar_fac, ftol_rel = control$ftol_rel)

  .create_PF_clouds(out)
}

.create_PF_clouds <- function(
  forward_clouds, backward_clouds = list(), smoothed_clouds = list(),
  transition_likelihoods = list())
  structure(list(
    forward_clouds = forward_clouds, backward_clouds = backward_clouds,
    smoothed_clouds = smoothed_clouds,
    transition_likelihoods = transition_likelihoods), class = "PF_clouds")

get_Q_tilde <- function(x, n_vars)
  if(is.null(x)) diag(0, n_vars) else x

#' @importFrom graphics plot
.PF_EM <- function(
  n_fixed_terms_in_state_vec, X, fixed_terms, tstart, tstop, Q_0, Q, a_0, F.,
  R, risk_obj, n_max, n_threads, N_fw_n_bw, N_smooth, N_smooth_final, N_first,
  eps, nu, covar_fac,
  forward_backward_ESS_threshold = NULL, debug = 0, trace,
  method = "AUX_normal_approx_w_particles", seed = NULL, smoother, model,
  fixed_params, type, Q_tilde, est_a_0, G, J, K, theta, psi, phi, ftol_rel,
  averaging_start, fix_seed){
  cl <- match.call()
  n_vars <- nrow(X)
  fit_call <- cl
  fit_call[[1]] <- as.name("PF_smooth")

  fit_call[["Q_tilde"]] <- get_Q_tilde(Q_tilde, n_vars)

  is_restricted <- all(
    !is.null(G), !is.null(J), !is.null(K), !is.null(theta), !is.null(psi),
    !is.null(phi))
  if(is_restricted){
    fit_call[["F"]]   <- F. <- .get_F(G, theta)
    fit_call[["Q"]]   <- Q  <- .get_Q(J, K, psi, phi)$Q
    fit_call[["a_0"]] <- eval(fit_call[["a_0"]], parent.frame())
    G_tilde <- .get_cum_mat(nrow(F.), ncol(F.)) %*% G
    J_qr <- qr(J)

  } else {
    fit_call[["F"]]   <- eval(fit_call[["F."]] , parent.frame())
    fit_call[["a_0"]] <- eval(fit_call[["a_0"]], parent.frame())
    fit_call[["Q"]]   <- eval(fit_call[["Q"]]  , parent.frame())

  }

  fit_call[c("eps", "seed", "F.", "trace", "est_a_0", "G", "J", "K",
             "theta", "psi", "phi", "averaging_start", "fix_seed")] <- NULL

  # print Q and F structure
  if(trace > 0 && is_restricted){
    tmp <- list(
      F = .get_F(G, seq_along(theta)), R = R,
      Q = .get_Q(J, K, seq_along(psi), seq_along(phi))$Q)
    tmp <- .set_PF_names(tmp, rng_names = row.names(X), fixed_names = NULL)

    # start with F
    if(all(rowSums(G) < 2, G %in% c(0, 1))){
      tmp$F <- structure(
        sapply(tmp$F, sprintf, fmt = "t%d"), dimnames = dimnames(tmp$F),
        dim = dim(tmp$F))
      tmp$F[tmp$F == "t0"] <- NA_character_
      cat(sQuote("F"), "matrix is of the following form\n")
      print(tmp$F, quote = FALSE, na.print = "")
      cat("\n")

    }

    if(all(rowSums(K) < 2, rowSums(J) < 2, c(J, K) %in% c(0, 1))){
      tp <- drop(J %*% seq_along(psi))
      diag(tmp$Q) <- sprintf("psi%d", tp)
      tp <- drop(K %*% seq_along(phi))
      tmp$Q[lower.tri(tmp$Q)] <- sprintf("phi%d", tp)
      tmp$Q[tmp$Q == "phi0"] <- tmp$Q[upper.tri(tmp$Q)] <- NA_character_
      cat(sQuote("Q"), "matrix is of the following form\n")
      print(tmp$Q, quote = FALSE, na.print = "")
      cat("\n")
    }

  }

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

  # setup matrices and arrays to save the estimates from each iteration
  a_0_it          <- matrix(NA_real_, n_max, length(a_0))
  fixed_params_it <- matrix(NA_real_, n_max, length(fixed_params))
  F_it <- array(NA_real_, c(nrow(F.), ncol(F.), n_max))
  Q_it <- array(NA_real_, c(nrow(Q) , ncol(Q) , n_max))

  if(type == "VAR")
    Q_0 <- fit_call$Q_0 <- get_Q_0(Qmat = fit_call$Q, Fmat = fit_call$F)

  log_likes <- rep(NA_real_, n_max)
  log_like <- log_like_max <- -Inf

  if(trace > 0){
    print_covmat <- function(X, lab){
      cormat <- cov2cor(X)
      cormat[upper.tri(cormat, diag = TRUE)] <- NA_real_
      cormat <- cormat[-1, -nrow(cormat)]
      cat("Square root of diagonal elements of", lab, "are:\n")
      print(sqrt(diag(X)), na.print = "")
      cat("Lower triangle of correlation matrix is\n")
      print(cormat, na.print = "")
    }
  }

  for(i in 1:n_max){
    if(trace > 0){
      if(i != 1)
        cat("\n")
      cat("#######################\nStarting EM iteration", i, "\n")

      cat("a_0 is:\n")
      print(fit_call$a_0)

      if(length(fixed_params) > 0){
        cat("Fixed parameters are:\n")
        print(fixed_params)
      }
      if(type == "VAR"){
        cat("F is:\n")
        print(fit_call$F)
      }

      print_covmat(fit_call$Q, "Q")

      if(type == "VAR")
        print_covmat(fit_call$Q_0, "Q_0")
    }
    log_like_old <- log_like
    log_like_max <- max(log_like, log_like_max)

    #####
    # find clouds
    if(fix_seed)
      assign(".Random.seed", seed, envir = .GlobalEnv)
    clouds <- eval(fit_call, envir = parent.frame())

    if(trace > 0){
      cat("Plotting state vector mean and quantiles for iteration", i,
          "(dashed lines are means from forward and backward particle",
          "clouds)\n")
      plot(clouds, main = paste0("EM iteration ", i), qlvls = c(.95, .05),
           qtype = "lines")
      plot(clouds, type = "forward_clouds", add = TRUE, qlvls = c(), lty = 2)
      plot(clouds, type = "backward_clouds", add = TRUE, qlvls = c(), lty = 3)

      cat("Effective sample sizes are:\n")
      effective_sample_size <- PF_effective_sample_size(clouds)
      if(length(effective_sample_size[[1L]] > 50)){
        print(sapply(effective_sample_size, function(x) c(
          `Mean effective sample size` = mean(x),
          `Sd of effective sample size` = sd(x),
          Min = min(x), Max = max(x))))

      } else
        print(effective_sample_size)
    }

    #####
    # update parameters in state equation
    if(trace > 0)
      cat("Updating parameters in state model...\n")
    a_0_old <- fit_call$a_0
    Q_old <- fit_call$Q
    F_old <- fit_call$F
    fixed_params_old <- fixed_params

    if(type == "RW"){
      sum_stats <- compute_PF_summary_stats(
        clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
        debug = trace > 2, F = F.)
      if(est_a_0)
        a_0 <- drop(sum_stats[[1]]$E_xs)
      Q <- Reduce(
        "+", lapply(sum_stats, "[[", "E_x_less_x_less_one_outers"))
      Q <- Q / length(sum_stats)

      fit_call$a_0 <- a_0
      fit_call$Q <- Q

    } else if (type == "VAR") {
      if(is_restricted){
        if(trace > 0)
          cat("Running first conditional maximization step\n")
        new_params <- PF_est_params_dens(
          clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
          debug = trace > 1, only_QR = TRUE)
        # see https://math.stackexchange.com/a/2875857/253239
        QR_R <- new_params$QR_R
        QR_F <- new_params$QR_F
        t1 <- crossprod(
          G_tilde, as.vector(crossprod(QR_R, t(solve(Q, t(QR_F))))))
        t2 <- crossprod(
          G_tilde, kronecker(solve(Q), crossprod(QR_R)) %*% G_tilde)
        theta <- drop(solve(t2, t1))
        # TODO: need to change for higher order models
        fit_call$F <- F. <- .get_F(G, theta)

        if(trace > 0)
          cat("Running second conditional maximization step\n")
        sum_stats <- compute_PF_summary_stats(
          clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
          debug = trace > 2, F = F., do_use_F = TRUE, do_compute_E_x = FALSE)
        Z <- Reduce(
          "+", lapply(sum_stats, "[[", "E_x_less_x_less_one_outers"))

        #####
        # see https://stats.stackexchange.com/q/362062/81865

        # assign log-likelihood function
        idx <- 1:length(psi)
        nobs <- length(sum_stats)
        ll <- function(par){
          psi <- par[ idx]
          phi <- par[-idx]

          Q <- .get_Q(J, K, psi, phi)$Q
          Q_qr <- qr(Q)
          deter <- determinant(Q, logarithm = TRUE)
          if(deter$sign < 0 || Q_qr$rank < ncol(Q))
            return(NA_real_)

          -(nobs * deter$modulus + sum(diag(solve(Q_qr, Z)))) / 2
        }

        # assign gradient function
        gr <- function(par){
          psi <- par[ idx]
          phi <- par[-idx]

          tmp <- .get_Q(J, K, psi, phi)
          Q <- tmp$Q
          C <- tmp$C
          V <- tmp$V

          # TODO: Computations could be done a lot smarter...
          Q_qr <- qr(Q)
          if(Q_qr$rank < ncol(Q) || det(Q) <= 0)
            return(NA_real_)
          fac <- solve(Q_qr, Z) - diag(nobs, ncol(Q))
          fac <- solve(Q_qr, t(fac)) / 2

          d_V <-
            diag(fac %*% V %*% C + C %*% V %*% fac) %*%
            J %*% diag(exp(psi), length(psi))

          d_C <- tcrossprod(V, V %*% fac)
          exp_phi <- exp(phi)
          d_C <- as.vector(d_C)[lower.tri(d_C, diag = FALSE)] %*% K %*%
            diag((4 * exp_phi / (1 + exp_phi)^2), length(phi))

          c(drop(d_V), drop(d_C))
        }

        # scale with full covariance matrix log-likelihood
        Q_full <- Z / nobs
        deter <- determinant(Q_full, logarithm = TRUE)
        ll_full <- -(nobs * deter$modulus + nobs * ncol(Q_full)) / 2
        ctrl <- list(fnscale = -abs(ll_full),
                     reltol = .Machine$double.eps^(3/4))

        out_optim <-
          optim(c(psi, phi), ll, gr, control = ctrl, method = "BFGS")
        if(out_optim$convergence != 0)
          stop(sQuote("optim"), " failed to converge with code ",
               out_optim$convergence)
        psi <- out_optim$par[ idx]
        phi <- out_optim$par[-idx]

        fit_call$Q <- Q <- .get_Q(J, K, psi, phi)$Q

      } else {
        new_params <- PF_est_params_dens(
          clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
          debug = trace > 1, only_QR = FALSE)
        # TODO: need to change for higher order models
        fit_call$F <- F. <- new_params$R_top_F
        fit_call$Q <- Q <- new_params$Q

      }


    } else
      stop(sQuote("type"), " not implemented")

    #####
    # Update fixed effects
    has_fixed_params <- length(fixed_params) > 0
    if(has_fixed_params){
      if(trace > 0)
        cat("Updating fixed effects...\n")

      fit_call$fixed_params <- fixed_params <- .PF_update_fixed(
        clouds = clouds$smoothed_clouds, risk_obj = risk_obj, model = model,
        R = R, X = X, fixed_terms = fixed_terms, fixed_params = fixed_params,
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

    a_0_it[i, ] <- a_0
    if(has_fixed_params)
      fixed_params_it[i, ] <- fixed_params
    F_it[, ,i] <- F.
    Q_it[, ,i] <- Q

    # average if requested
    if(averaging_start > 0 && i > averaging_start){
      avg_idx <- averaging_start:i
      # TODO: could be done much faster
      fit_call$a_0 <- a_0 <- colMeans(a_0_it[avg_idx, , drop = FALSE])
      if(has_fixed_params)
        fit_call$fixed_params <- fixed_params <- colMeans(
          fixed_params_it[avg_idx, , drop = FALSE])

      fit_call$F <- F. <- apply(F_it[, , avg_idx, drop = FALSE], 1:2, mean)
      fit_call$Q <- Q  <- apply(Q_it[, , avg_idx, drop = FALSE], 1:2, mean)

    }

    if(type == "VAR")
      Q_0 <- fit_call$Q_0 <- get_Q_0(Qmat = fit_call$Q, Fmat = fit_call$F)

    # compute norms
    Q_relative_norm <- norm(Q_old - Q) / (norm(Q_old) + 1e-8)
    a_0_relative_norm <- norm(t(a_0 - a_0_old)) / (norm(t(a_0_old)) + 1e-8)
    F_norm <- norm(F_old - F.) / (norm(F_old) + 1e-8)
    if(has_fixed_params)
      fixed_params_norm <- norm(t(fixed_params - fixed_params_old)) /
      (norm(t(fixed_params)) + 1e-8)

    if(trace > 0){
      msg <- "The relative norm of the change in"
      if(type == "VAR")
        msg <- paste(msg, "F,")

      if(has_fixed_params)
        msg <- paste(msg, "fixed parameters,")

      msg <- paste(msg, "a_0 and Q are")
      if(type == "VAR")
        msg <- paste0(msg, " ", sprintf("%.5f", F_norm), ",")

      if(has_fixed_params)
        msg <- paste0(msg, " ", sprintf("%.5f", fixed_params_norm), ",")

      msg <- paste(
        msg, sprintf("%.5f", a_0_relative_norm), "and",
        sprintf("%.5f", Q_relative_norm), "at iteration", i)
      cat(msg, "\n")
    }

    if(has_converged <-
       Q_relative_norm < eps &&
       a_0_relative_norm < eps &&
       (type != "VAR" || F_norm < eps) &&
       (!has_fixed_params || fixed_params_norm < eps))
      break
  }

  if(!has_converged)
    warning("Method did not converge.")

  if(!exists("effective_sample_size", envir = environment()))
    effective_sample_size <- PF_effective_sample_size(clouds)

  out <- structure(list(
    call = cl, clouds = clouds, a_0 = a_0, fixed_effects = fixed_params, Q = Q,
    F = fit_call$F, R = R, EM_ests = list(
      a_0             = a_0_it         [1:i, , drop = FALSE],
      fixed_effects   = fixed_params_it[1:i, , drop = FALSE],
      F = F_it[, , 1:i, drop = FALSE],
      Q = Q_it[, , 1:i, drop = FALSE]),
    log_likes = log_likes[1:i], n_iter = i,
    effective_sample_size = effective_sample_size, seed = seed),
    class = "PF_EM")

  if(is_restricted)
    out[c("G", "J", "K", "theta", "psi", "phi")] <- list(
      G = G, J = J, K = K, theta = theta, psi = psi, phi = phi)

  out
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
#' @param N_smooth_final number of particles to sample with replacement from
#' the smoothed particle cloud with \code{N_smooth} particles using the
#' particles' weights. This causes additional sampling error but decreases the
#' computation time in the M-step.
#' @param nu integer with degrees of freedom to use in the (multivariate)
#' t-distribution used as the proposal distribution. A (multivariate) normal
#' distribution is used if it is zero.
#' @param covar_fac factor to scale the covariance matrix with. Ignored if
#' the values is less than or equal to zero.
#' @param ftol_rel relative convergence tolerance of the mode objective in mode
#' approximation.
#' @param averaging_start index to start averaging. Values less then or equal
#' to zero yields no averaging.
#' @param fix_seed \code{TRUE} if the same seed should be used. E.g., in
#' \code{\link{PF_EM}} the same seed will be used in each iteration of the
#' E-step of the MCEM algorithm.
#'
#' @details
#' The \code{method} argument can take the following values
#'
#' \itemize{
#' \item \code{bootstrap_filter} for a bootstrap filter.
#' \item \code{PF_normal_approx_w_cloud_mean} for a particle filter where a
#' Gaussian approximation is used using a Taylor
#' approximation made at the mean for the current particle given the mean of the
#' parent particles  and/or mean of the child particles.
#' \item \code{AUX_normal_approx_w_cloud_mean} for an auxiliary particle filter
#' version of \code{PF_normal_approx_w_cloud_mean}.
#' \item \code{PF_normal_approx_w_particles} for a filter similar to
#' \code{PF_normal_approx_w_cloud_mean} and differs by making a Taylor
#' approximation at a mean given each sampled parent and/or child particle.
#' \item \code{AUX_normal_approx_w_particles} for an auxiliary particle filter
#' version of \code{PF_normal_approx_w_particles}.
#' }
#'
#' The \code{smoother} argument can take the following values
#' \itemize{
#' \item \code{Fearnhead_O_N} for the smoother in Fearnhead, Wyncoll, and Tawn
#' (2010).
#' \item \code{Brier_O_N_square} for the smoother in Briers, Doucet, and
#' Maskell (2010).
#' }
#'
#' @references
#' Gordon, N. J., Salmond, D. J., and Smith, A. F. (1993) Novel approach
#' to nonlinear/non-Gaussian Bayesian state estimation.
#' \emph{In IEE Proceedings F (Radar and Signal Processing)},
#' (Vol. 140, No. 2, pp. 107-113). IET Digital Library.
#'
#' Pitt, M. K., and Shephard, N. (1999) Filtering via simulation: Auxiliary
#' particle filters. \emph{Journal of the American statistical association},
#' \strong{94(446)}, 590-599.
#'
#' Fearnhead, P., Wyncoll, D., and Tawn, J. (2010) A sequential smoothing
#' algorithm with linear computational cost. \emph{Biometrika}, \strong{97(2)},
#' 447-464.
#'
#' Briers, M., Doucet, A., and Maskell, S. (2010) Smoothing algorithms for
#' state-space models.
#' \emph{Annals of the Institute of Statistical Mathematics}, \strong{62(1)},
#' 61.
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
  Q_tilde = NULL, est_a_0 = TRUE, N_smooth_final = N_smooth, nu = 0L,
  covar_fac = -1, ftol_rel = 1e-8, averaging_start = -1L, fix_seed = TRUE){
  control <- list(
    N_fw_n_bw = N_fw_n_bw, N_smooth = N_smooth, N_first = N_first, eps = eps,
    forward_backward_ESS_threshold = forward_backward_ESS_threshold,
    method = method, n_max = n_max, n_threads = n_threads, smoother = smoother,
    Q_tilde = Q_tilde, est_a_0 = est_a_0, N_smooth_final = N_smooth_final,
    nu = nu, covar_fac = covar_fac, ftol_rel = ftol_rel,
    averaging_start = averaging_start, fix_seed = fix_seed)

  stopifnot(
    length(method) == 1L, method %in% c(
      "bootstrap_filter", "PF_normal_approx_w_cloud_mean",
      "AUX_normal_approx_w_cloud_mean", "PF_normal_approx_w_particles",
      "AUX_normal_approx_w_particles"),
    length(smoother) == 1L, smoother %in% c(
      "Fearnhead_O_N", "Brier_O_N_square"))

  check_n_particles_expr <- function(N_xyz)
    eval(bquote({
      if(is.null(control[[.(N_xyz)]]))
        stop("Please supply the number of particle in ", sQuote(paste0(
          "PF_control(", .(N_xyz), ")")))
    }), parent.frame())

  check_n_particles_expr("N_first")
  check_n_particles_expr("N_fw_n_bw")
  check_n_particles_expr("N_smooth")
  stopifnot(
    typeof(N_smooth_final) %in% c("double", "integer"),
    length(N_smooth_final) == 1L, as.integer(N_smooth_final) == N_smooth_final,
    N_smooth_final <= N_smooth,
    typeof(nu) %in% c("double", "integer"), length(nu) == 1L, nu >= 0L,
    as.integer(nu) == nu, is.numeric(covar_fac), is.numeric(ftol_rel),
    ftol_rel > 0,
    is.integer(averaging_start),
    is.logical(fix_seed), length(fix_seed) == 1L)

  return(control)
}


.get_PF_static_args <- function(
  formula, data, by, max_T = NULL, id, trace, model, order, fixed = NULL,
  random = NULL){
  # get design matrix and risk set
  tmp <- get_design_matrix_and_risk_obj(
    formula = formula, data = data, by = by,
    max_T = if(is.null(max_T)) NULL else max_T, verbose = trace > 0,
    is_for_discrete_model = model %in% c("logit", "cloglog"), id = id,
    fixed = fixed, random = random)

  if(trace > 0)
    report_pre_liminary_stats_before_EM(
      risk_set = tmp$risk_set, Y = tmp$X_Y$Y)

  # tranpose due to column-major storage and we want to look up individuals
  tmp$X_Y$X           <- t(tmp$X_Y$X)
  tmp$X_Y$fixed_terms <- t(tmp$X_Y$fixed_terms)

  with(tmp, list(
    n_fixed_terms_in_state_vec = 0, X = X_Y$X, fixed_terms = X_Y$fixed_terms,
    tstart = X_Y$Y[, 1], tstop = X_Y$Y[, 2], risk_obj = risk_set,
    debug = max(0, trace - 1), model = model, terms = X_Y$terms,
    has_fixed_intercept = X_Y$has_fixed_intercept,
    xlev = X_Y$xlev))
}

get_family_arg <- function(model)
  switch(
    model, logit = "binomial", cloglog = "cloglog", exponential = "poisson")

.PF_update_fixed <- function(
  clouds, risk_obj, R, X, fixed_terms, fixed_params, model, nthreads,
  tstart, tstop, debug){
  if(!model %in% c("logit", "cloglog", "exponential"))
    stop(sQuote(model), " is not implemented with fixed effects")

  family_arg <- get_family_arg(model)

  R_top <- t(R)
  out <- pf_fixed_effect_get_QR(
    clouds = clouds, risk_obj = risk_obj, ran_vars = X,
    fixed_terms = fixed_terms, R_top = R_top, tstart = tstart,
    tstop = tstop, fixed_params = fixed_params, family = family_arg,
    max_threads = nthreads, debug = debug)

  Xty <- Reduce("+", lapply(out, "[[", "XtWY"))
  qr_o <- qr(do.call(rbind, lapply(out, "[[", "Rmat")))
  qr_R <- qr.R(qr_o)
  drop(solve(qr_R, solve(t(qr_R), Xty)))
}

.get_Q <- function(J, K, psi, phi){
  V <- diag(exp(drop(J %*% psi)))
  C <- diag(1, ncol(V))
  C[lower.tri(C)] <- 2/(1 + exp(-drop(K %*% phi))) - 1
  C[upper.tri(C)] <- t(C)[upper.tri(C)]
  list(Q = V %*% C %*% V, V = V, C = C)
}

.get_F <- function(G, theta){
  q <- as.integer(sqrt(nrow(G)) + 1e-8)
  matrix(as.vector(G %*% theta), q, q)
}

.get_cum_mat <- function(r, c){
  Ir <- diag(1, r)
  Ic <- diag(1, c)
  H <- list()
  for (i in 1:r) {
    H[[i]] <- list()
    for (j in 1:c) {
      H[[i]][[j]] <- Ir[i, ] %o% Ic[j, ]
    }
  }
  p <- r * c
  K <- matrix(0, nrow = p, ncol = p)
  for (i in 1:r) {
    for (j in 1:c) {
      Hij <- H[[i]][[j]]
      K <- K + (Hij %x% t(Hij))
    }
  }
  K
}

# see https://gist.github.com/boennecd/09ab5b0baae4738089530ae37bc9812e
.get_dup_mat <- function(n){
  if(n == 1L)
    return(as.matrix(1.))

  stopifnot(is.integer(n), n > 1L)
  nq <- n * (n + 1L) / 2L

  o <- matrix(0L, n, n)
  o[lower.tri(o, diag = TRUE)] <- 1:nq
  o[upper.tri(o)] <- t(o)[upper.tri(o)]
  o <- c(o)

  nn <- n * n
  out <- matrix(0L, nn, nq)
  for(i in 1:nn)
    out[i, o[i]] <- 1L

  out
}


#' @title Compute Time-Invariant Covariance Matrix
#' @description
#' Computes the invariant covariance matrix for a vector autoregression model.
#'
#' @param Qmat covariance matrix in transition density.
#' @param Fmat coefficients in transition density.
#'
#' @return
#' The invariant covariance matrix.
#'
#' @examples
#' Fmat <- matrix(c(.8, .4, .1, .5), 2, 2)
#' Qmat <- matrix(c( 1, .5, .5,  2), 2)
#'
#' x1 <- get_Q_0(Qmat = Qmat, Fmat = Fmat)
#' x2 <- Qmat
#' for(i in 1:101)
#'   x2 <- tcrossprod(Fmat %*% x2, Fmat) + Qmat
#' stopifnot(isTRUE(all.equal(x1, x2)))
#'
#'
#' @export
get_Q_0 <- function(Qmat, Fmat){
  eg  <- eigen(Fmat)
  las <- eg$values
  if(any(Mod(las) >= 1))
    stop("Divergent series")
  U   <- eg$vectors
  T. <- solve(U, t(solve(U, Qmat)))
  Z   <- T. / (1 - tcrossprod(las))
  out <- tcrossprod(U %*% Z, U)
  if(is.complex(out)){
    if(all(abs(Im(out)) < .Machine$double.eps^(3/4)))
      return(Re(out))

    stop("Q_0 has imaginary part")
  }

  out
}

#' @name get_cloud_means
#' @title Compute Mean Estimates from Particle Cloud
#' @description
#' Computes the estimated means from a particle cloud.
#'
#' @param object object with class \code{PF_EM} or \code{PF_clouds}.
#' @param cov_index integer vector with indices of the random effect to
#' include.
#' @param type character with the type of cloud to compute means for.
#' @param ... named arguments to pass to the \code{PF_clouds} method.
#'
#' @return
#' A matrix which rows are time indices and columns are random effect indices.
#'
#' @export
get_cloud_means <- function(object, ...){
  UseMethod("get_cloud_means")
}

#' @rdname get_cloud_means
#' @export
get_cloud_means.PF_EM <- function(object, ...){
  cl <- match.call()
  cl[[1]] <- quote(get_cloud_means)
  cl$object <- bquote(.(substitute(object))$clouds)
  eval(cl, parent.frame())
}

#' @rdname get_cloud_means
#' @export
get_cloud_means.PF_clouds <- function(
  object, cov_index = NULL,
  type = c("smoothed_clouds", "forward_clouds", "backward_clouds"),
  ...)
{
  type <- type[1L]
  stopifnot(
    type %in% c("smoothed_clouds", "forward_clouds", "backward_clouds"))

  cl <- object[[type]]

  if(is.null(cov_index))
    cov_index <- seq_len(dim(cl[[1L]]$states)[1])

  do.call(rbind, sapply(cl, function(row){
    colSums(t(row$states[cov_index, , drop = FALSE]) * drop(row$weights))
  }, simplify = FALSE))
}





#' @name get_cloud_quantiles
#' @title Compute Quantile Estimates from Particle Cloud
#' @description
#' Computes the estimated quantiles from a particle cloud.
#'
#' @inheritParams get_cloud_means
#' @param type character with the type of cloud to compute quantiles for.
#' @param qlvls numeric vector with values in \eqn{[0,1]} with the quantiles to
#' compute.
#'
#' @return
#' A 3 dimensional array where the first dimension is the quantiles, the second
#' dimension is the random effect, and the third dimension is the time.
#'
#' @export
get_cloud_quantiles <- function(object, ...){
  UseMethod("get_cloud_quantiles")
}

#' @rdname get_cloud_quantiles
#' @export
get_cloud_quantiles.PF_EM <- function(object, ...){
  cl <- match.call()
  cl[[1]] <- quote(get_cloud_quantiles)
  cl$object <- bquote(.(substitute(object))$clouds)
  eval(cl, parent.frame())
}

#' @rdname get_cloud_quantiles
#' @export
get_cloud_quantiles.PF_clouds <- function(
  object, cov_index = NULL, qlvls = c(.05, .5, .95),
  type = c("smoothed_clouds", "forward_clouds", "backward_clouds"),
  ...)
{
  type <- type[1L]
  stopifnot(
    type %in% c("smoothed_clouds", "forward_clouds", "backward_clouds"))
  stopifnot(length(qlvls) > 0, all(0 <= qlvls, qlvls <= 1))

  cl <- object[[type]]

  if(is.null(cov_index))
    cov_index <- seq_len(dim(cl[[1L]]$states)[1])

  qs <- lapply(cl, function(row){
    out <- apply(row$states[cov_index, , drop = FALSE], 1, function(x){
      ord <- order(x)
      wg_cumsum <- cumsum(row$weights[ord])
      idx <- ord[sapply(qlvls, function(q) {
        is_lower <- wg_cumsum < q
        if(!any(is_lower))
          return(NA_integer_)
        max(which(wg_cumsum < q))
      })]
      x[idx]
    })

    if(is.null(dim(out)))
      out <- matrix(out, ncol = length(out))

    out
  })

  if(identical(dim(qs[[1]]), c(1L, 1L)))
    return(array(simplify2array(qs), dim = c(1L, 1L, length(qs))))

  simplify2array(qs)
}


#' @title Approximate Observed Information Matrix and Score Vector
#' @description
#' Returns a list of functions to approximate the observed information matrix
#' and score vector.
#'
#' @param object object of class \code{\link{PF_EM}}.
#' @param debug \code{TRUE} if debug information should be printed to the
#' console.
#' @param use_O_n_sq \code{TRUE} if the method from Poyiadjis et al. (2011)
#' should be used.
#'
#' @details
#' The score vector and observed information matrix are computed
#' with the (forward)
#' particle filter. This comes at an \eqn{O(d^2)} variance where \eqn{d}
#' is the number of periods. Thus, the approximation may be poor for long
#' series. The score vector can be used to perform stochastic gradient
#' descent.
#'
#' If \code{use_O_n_sq} is \code{TRUE} then the method in Poyiadjis et al. (2011)
#' is used. This may only have a variance which is linear in the number of
#' time periods. However, the present implementation is \eqn{O(N^2)} where
#' \eqn{N} is the number of particles. The method uses a particle filter as
#' in Section 3.1
#' of Lin et al. (2005). There is no need to call
#' \code{run_particle_filter} unless one wants a new approximation of the
#' log-likelihood as a separate filter is run with \code{get_get_score_n_hess}
#' when \code{use_O_n_sq} is \code{TRUE}.
#'
#' @section Warning:
#' The function is still under development so the output and API may change.
#'
#' @seealso
#' See the examples at https://github.com/boennecd/dynamichazard/tree/master/examples.
#'
#' @references
#' Cappe, O. and Moulines, E. (2005) Recursive Computation of the Score and
#' Observed Information Matrix in Hidden Markov Models.
#' \emph{IEEE/SP 13th Workshop on Statistical Signal Processing}.
#'
#' Cappe, O., Moulines, E. and Ryden, T. (2005) Inference in Hidden Markov
#' Models (Springer Series in Statistics). Springer-Verlag.
#'
#' Doucet, A., and Tadić, V. B. (2003) Parameter Estimation in General
#' State-Space Models Using Particle Methods.
#' \emph{Annals of the Institute of Statistical Mathematics}, \strong{55(2)},
#' 409–422.
#'
#' Lin, M. T., Zhang, J. L., Cheng, Q. and Chen, R. (2005) Independent
#' Particle Filters. \emph{Journal of the American Statistical Association},
#' \strong{100(472)}, 1412-1421.
#'
#' Poyiadjis, G., Doucet, A. and Singh, S. S. (2011) Particle Approximations of
#' the Score and Observed Information Matrix in State Space Models with
#' Application to Parameter Estimation. \emph{Biometrika}, \strong{98(1)},
#' 65--80.
#'
#' @return
#' A list with the following functions as elements
#' \item{run_particle_filter}{function to run particle filter as with
#' \code{\link{PF_forward_filter}}.}
#' \item{set_parameters}{function to set the parameters in the model.
#' The first argument is a vectorized version of \eqn{F} matrix and \eqn{Q}
#' matrix. The second argument is the fixed effect coefficients.}
#' \item{set_n_particles}{sets the number of particles to use in
#' \code{run_particle_filter} and \code{get_get_score_n_hess} when
#' \code{use_O_n_sq} is \code{TRUE}.}
#' \item{get_get_score_n_hess}{approximate the observed information
#' matrix and score vector. The argument toggles whether or not to approximate
#' the observed information matrix. The last particle cloud
#' from \code{run_particle_filter} is used when \code{use_O_n_sq} is
#' \code{FALSE}.}
#'
#' @examples
#' \dontrun{
#' library(dynamichazard)
#' .lung <- lung[!is.na(lung$ph.ecog), ]
#' # standardize
#' .lung$age <- scale(.lung$age)
#'
#' # fit model
#' set.seed(43588155)
#' pf_fit <- PF_EM(
#'   fixed = Surv(time, status == 2) ~ ph.ecog + age,
#'   random = ~ 1, model = "exponential",
#'   data = .lung, by = 50, id = 1:nrow(.lung),
#'   Q_0 = as.matrix(1), Q = as.matrix(.5^2), type = "VAR",
#'   max_T = 800, Fmat = as.matrix(.5),
#'   control = PF_control(
#'     N_fw_n_bw = 250, N_first = 2000, N_smooth = 500, covar_fac = 1.1,
#'     nu = 6, n_max = 1000L, eps = 1e-4, averaging_start = 200L,
#'     n_threads = max(parallel::detectCores(logical = FALSE), 1)))
#'
#' # compute score and observed information matrix
#' comp_obj <- PF_get_score_n_hess(pf_fit)
#' comp_obj$set_n_particles(N_fw = 10000L, N_first = 10000L)
#' comp_obj$run_particle_filter()
#' (o1 <- comp_obj$get_get_score_n_hess())
#'
#' # O(N^2) method with lower variance as a function of time
#' comp_obj <- PF_get_score_n_hess(pf_fit, use_O_n_sq = TRUE)
#' comp_obj$set_n_particles(N_fw = 2500L, N_first = 2500L)
#' (o2 <- comp_obj$get_get_score_n_hess())
#'
#' # approximations may have large variance
#' o3 <- replicate(10L, {
#'   runif(1)
#'   pf_fit$seed <- .Random.seed
#'   comp_obj <- PF_get_score_n_hess(pf_fit)
#'   comp_obj$set_n_particles(N_fw = 10000L, N_first = 10000L)
#'   comp_obj$run_particle_filter()
#'   comp_obj$get_get_score_n_hess()
#' }, simplify = FALSE)
#' sapply(o3, function(x) x$score)
#' sapply(o3, function(x) sqrt(diag(solve(x$obs_info))))
#' }
#' @export
PF_get_score_n_hess <- function(object, debug = FALSE, use_O_n_sq = FALSE){
  stopifnot(inherits(object, "PF_EM"))

  #####
  # get design matrix, etc. again
  org_cl <- object$call
  ma <- match(
    c("formula", "data", "by", "max_T", "id", "trace", "model", "order",
      "fixed", "random"), names(org_cl), nomatch = 0L)
  sta_arg_call <- org_cl[c(1L, ma)]
  if("data" %in% names(sta_arg_call))
    cat("Using",
        sQuote(paste0(deparse(sta_arg_call$data), collapse = "\n")),
        "as the", sQuote("data"), "argument\n")
  if(!"id" %in% names(sta_arg_call)){
    # TODO: check that this works
    if(!"data" %in% names(sta_arg_call))
      stop(sQuote("data"), " is needed when id is not used in original call")

    sta_arg_call$id <- bquote(seq_len(nrow(.(sta_arg_call$data))))
  }

  # add defaults where needed
  def_args <- c("model", "fixed", "random", "order")
  if(any(is_missing <- !def_args %in% names(sta_arg_call)))
    sta_arg_call[def_args[is_missing]] <- formals(PF_EM)[def_args[is_missing]]
  sta_arg_call$trace <- FALSE
  sta_arg_call[[1L]] <- quote(list)

  static_args <-
    do.call(.get_PF_static_args, eval(sta_arg_call, parent.frame()))

  # find model argument
  model <- if(is.null(org_cl$model)) formals(PF_EM)$model else org_cl$model
  family_arg <- get_family_arg(model)

  # handle additional arguments
  ctrl <- object$control
  seed <- object$seed
  type <- if(is.null(object$call$type))
    formals(PF_EM)$type else object$call$type

  fixed_effects <- object$fixed_effects
  Q <- object$Q
  n_q <- length(Q)
  Fmat <- object$F
  n_f <- length(Fmat)
  R <- object$R
  a_0 <- object$a_0
  Q_0 <- if(type == "VAR") get_Q_0(Q, Fmat) else eval(object$call$Q_0)
  if(!is.matrix(Q_0))
    Q_0 <- as.matrix(Q_0)

  N_fw <- ctrl$N_fw_n_bw
  N_first <- ctrl$N_first
  Q_tilde <- get_Q_tilde(ctrl$Q_tilde, ncol(Q))

  #####
  # define functions to return
  fw_cloud <- object$clouds$forward_clouds
  # runs particle filter and returns the particle clouds
  run_particle_filter <- function(){
    if(ctrl$fix_seed)
      assign(".Random.seed", seed, envir = .GlobalEnv)
    fw_cloud <<- particle_filter(
      fixed_params = fixed_effects, type = type, n_fixed_terms_in_state_vec =
        static_args$n_fixed_terms_in_state_vec, X = static_args$X,
      fixed_terms = static_args$fixed_terms, tstart = static_args$tstart,
      tstop = static_args$tstop, risk_obj = static_args$risk_obj,
      debug = debug, model = static_args$model, Q = Q, Q_0 = Q_0,
      F = Fmat, R = R, is_forward = TRUE, a_0 = a_0, N_fw_n_bw = N_fw,
      N_first = N_first, nu = if(is.null(ctrl$nu)) 0L else ctrl$nu,
      forward_backward_ESS_threshold = ctrl$forward_backward_ESS_threshold,
      method = ctrl$method, n_threads = ctrl$n_threads, Q_tilde = Q_tilde,
      covar_fac = ctrl$covar_fac, ftol_rel = ctrl$ftol_rel)

    .create_PF_clouds(fw_cloud)
  }

  # set the number of particles to use
  set_n_particles <- function(N_fw, N_first){
    if(!missing(N_fw))
      N_fw <<- N_fw
    if(!missing(N_first))
      N_first <<- N_first

    invisible()
  }

  # set the parameters in the model
  set_parameters <- function(state, obs){
    if(!missing(state)){
      Fmat[] <<- state[1:n_f]
      Q[] <<- state[n_f + 1:n_q]
    }
    if(!missing(obs))
      fixed_effects[] <<- obs

    if(type == "VAR")
      Q_0 <- get_Q_0(Q, Fmat)

    invisible(list(Fmat = Fmat, Q = Q, fixed_effects = fixed_effects,
                   Q_0 = Q_0))
  }

  # returns the score and potentially the negative Hessian estimates
  get_get_score_n_hess <- function(only_score = FALSE){
    if(ctrl$fix_seed)
      assign(".Random.seed", seed, envir = .GlobalEnv)
    cpp_res <- PF_get_score_n_hess_cpp(
      fw_cloud = fw_cloud, Q = Q, F = Fmat,
      risk_obj = static_args$risk_obj, ran_vars = static_args$X,
      fixed_terms = static_args$fixed_terms, tstart = static_args$tstart,
      tstop = static_args$tstop, fixed_params = fixed_effects,
      max_threads = ctrl$n_threads, family = family_arg,
      debug = debug, only_score = only_score, a_0 = a_0, R = R,
      Q_0 = Q_0, Q_tilde = Q_tilde, N_fw_n_bw = N_fw, N_first = N_first,
      nu = ctrl$nu, covar_fac = ctrl$covar_fac, ftol_rel = ctrl$ftol_rel,
      method = ctrl$method,
      forward_backward_ESS_threshold = ctrl$forward_backward_ESS_threshold,
      use_O_n_sq = use_O_n_sq)

    # compute observed information matrix and score vector. We have
    # to: multiply parts of them part it by a commutation, and parts by
    # a duplication matrix

    # first define a few lengths
    dfix <- length(fixed_effects)
    n_rng <- NCOL(Q)
    drng <- 2L * n_rng * n_rng
    org_dim <- dfix + drng

    # output dimension and matrix to multiply objects by
    out_dim <- dfix + n_rng * n_rng + (n_rng * (n_rng + 1L)) / 2L
    trans_mat <- matrix(0, out_dim, org_dim)

    ifix <- 1:dfix
    trans_mat[ifix, ifix]          <- diag(dfix)
    # TODO: avoid this by changing other computations in c++
    idF  <- dfix + 1:(n_rng * n_rng)
    trans_mat[idF, idF]            <- .get_cum_mat(n_rng, n_rng)
    idQ <- dfix + n_rng * n_rng + 1:((n_rng * (n_rng + 1L)) / 2L)
    trans_mat[idQ, -c(ifix, idF)] <- t(.get_dup_mat(n_rng))

    # transform objects
    score <- cpp_res$score <- drop(trans_mat %*% cpp_res$score)
    cpp_res$score_outer <- tcrossprod(
      trans_mat %*% cpp_res$score_outer, trans_mat)
    cpp_res$hess_terms <- tcrossprod(
      trans_mat %*% cpp_res$hess_terms, trans_mat)

    # compute output and return
    obs_info <- with(cpp_res,
                     tcrossprod(score) - score_outer - hess_terms)

    # set names
    fnames <- rownames(static_args$fixed_terms)
    rnames <- rownames(static_args$X)
    dnames <- c(fnames,
                outer(rnames, rnames, function(x, y) paste0("F:", x, ".", y)))
    tmp <- outer(rnames, rnames, function(x, y) paste0("Q:", x, ".", y))
    dnames <- c(dnames, tmp[lower.tri(tmp, diag = TRUE)])
    names(score) <- dnames
    dimnames(obs_info) <- list(dnames, dnames)

    list(score = score, obs_info = obs_info)
  }

  list(
    run_particle_filter = run_particle_filter,
    set_n_particles = set_n_particles,
    get_get_score_n_hess = get_get_score_n_hess,
    set_parameters = set_parameters)
}
