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
#' @param seed seed to set at the start of every EM iteration. See
#' \code{\link{set.seed}}.
#' @param type type of state model. Either \code{"RW"} for a [R]andom [W]alk or
#' "VAR" for [V]ector [A]uto[R]egression.
#' @param Fmat starting value for \eqn{F} when \code{type = "VAR"}. See
#' 'Details'.
#' @param fixed_effects starting values for fixed effects if any. See
#' \code{\link{ddFixed}}.
#' @param G,theta,J,K,psi,phi parameters for a restricted \code{type = "VAR"} model.
#' See the vignette mentioned in 'Details' and the examples linked to in
#' 'See Also'.
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
#' \eqn{\alpha_t}, are related to the output throught the linear predictors
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
#' The function is still under development so the ouput and API may change.
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
#'     Q_tilde = diag(.05^2, p),
#'     n_max = 100L,  # should maybe be larger
#'     smoother = "Fearnhead_O_N", eps = 1e-4,
#'     n_threads = 4L # depends on your cpu(s)
#'   ),
#'   trace = 1L)
#' plot(fit$log_likes) # log-likelihood approximation at each iterations
#'
#' # take more iterations with more particles
#' cl <- fit$call
#' ctrl <- cl[["control"]]
#' ctrl[c("N_fw_n_bw", "N_smooth", "N_first", "n_max")] <- list(
#'   400L, 1000L, 1000L, 25L)
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
  fixed_parems <- start_coefs$fixed_parems_start

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
    trace = trace, seed = seed, fixed_parems = fixed_parems,
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
    Q_tilde = control$Q_tilde, est_a_0 = control$est_a_0)

  out <- .set_PF_names(out, rng_names = row.names(static_args$X),
                       fixed_names = rownames(static_args$fixed_terms))

  out$call <- match.call()
  out
}

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

#' @describeIn PF_forward_filter Forward particle filter using the
#' estimates of an \code{\link{PF_EM}} call.
#' @export
PF_forward_filter.PF_EM <- function(x, N_fw, N_first, seed, ...){
  cl <- x$call
  cl <- cl[c(1, match(
    c("formula", "data", "model", "by", "max_T", "id", "control",
      formals(PF_control), "type", "trace", "Q_0"),
    names(cl), 0))]
  names(cl)[2] <- "x"
  xSym <- substitute(x)
  cl[c("seed", "Fmat", "a_0", "Q", "R", "fixed_effects")] <-
    lapply(
      c("seed",    "F", "a_0", "Q", "R", "fixed_effects"),
      function(z) substitute(y$z, list(y = xSym, z = as.symbol(z))))
  if(!missing(seed))
    cl[["seed"]] <- substitute(seed)

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
    N_first = N_first, nu = control$nu,
    forward_backward_ESS_threshold = control$forward_backward_ESS_threshold,
    method = control$method, n_threads = control$n_threads, Q_tilde = Q_tilde)

  structure(list(
    forward_clouds = out, backward_clouds = list(), smoothed_clouds = list(),
    transition_likelihoods = list()), class = "PF_clouds")
}

#' @importFrom graphics plot
.PF_EM <- function(
  n_fixed_terms_in_state_vec, X, fixed_terms, tstart, tstop, Q_0, Q, a_0, F.,
  R, risk_obj, n_max, n_threads, N_fw_n_bw, N_smooth, N_smooth_final, N_first,
  eps, nu,
  forward_backward_ESS_threshold = NULL, debug = 0, trace,
  method = "AUX_normal_approx_w_particles", seed = NULL, smoother, model,
  fixed_parems, type, Q_tilde, est_a_0, G, J, K, theta, psi, phi){
  cl <- match.call()
  n_vars <- nrow(X)
  fit_call <- cl
  fit_call[[1]] <- as.name("PF_smooth")

  if(is.null(Q_tilde))
    fit_call[["Q_tilde"]] <- diag(0, n_vars)

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
             "theta", "psi", "phi")] <- NULL

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
  fixed_parems_it <- matrix(NA_real_, n_max, length(fixed_parems))
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

      if(length(fixed_parems) > 0){
        cat("Fixed parameters are:\n")
        print(fixed_parems)
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
    assign(".Random.seed", seed, envir = .GlobalEnv)
    clouds <- eval(fit_call, envir = parent.frame())

    if(trace > 0){
      cat("Plotting state vector mean and quantiles for iteration", i,
          "(dashed lines are means from forward and backward particle",
          "clouds)\n")
      plot(clouds, main = paste0("EM iteration ", i))
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
    fixed_parems_old <- fixed_parems

    if(type == "RW"){
      sum_stats <- compute_PF_summary_stats(
        clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0, R = R,
        debug = trace > 2, F = F.)
      if(est_a_0)
        a_0 <- drop(sum_stats[[1]]$E_xs)
      Q <- Reduce(
        "+", lapply(sum_stats, "[[", "E_x_less_x_less_one_outers")[-1])
      Q <- Q / (length(sum_stats) - 1)

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
          "+", lapply(sum_stats, "[[", "E_x_less_x_less_one_outers")[-1])

        #####
        # see https://stats.stackexchange.com/q/362062/81865

        # assign log-likelihood function
        idx <- 1:length(psi)
        nobs <- length(sum_stats) - 1
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
        ctrl <- list(fnscale = if(ll_full < 0) ll_full else -ll_full,
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

    if(type == "VAR")
      Q_0 <- fit_call$Q_0 <- get_Q_0(Qmat = fit_call$Q, Fmat = fit_call$F)

    #####
    # Update fixed effects
    has_fixed_params <- length(fixed_parems) > 0
    if(has_fixed_params){
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
    if(has_fixed_params)
      fixed_params_norm <- norm(t(fixed_parems - fixed_parems_old)) /
        (norm(t(fixed_parems)) + 1e-8)

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

    a_0_it[i, ] <- a_0
    if(has_fixed_params)
      fixed_parems_it[i, ] <- fixed_parems
    F_it[, ,i] <- F.
    Q_it[, ,i] <- Q

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
    call = cl, clouds = clouds, a_0 = a_0, fixed_effects = fixed_parems, Q = Q,
    F = fit_call$F, R = R, EM_ests = list(
      a_0             = a_0_it         [1:i, , drop = FALSE],
      fixed_effects   = fixed_parems_it[1:i, , drop = FALSE],
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
  Q_tilde = NULL, est_a_0 = TRUE, N_smooth_final = N_smooth, nu = 0L){
  control <- list(
    N_fw_n_bw = N_fw_n_bw, N_smooth = N_smooth, N_first = N_first, eps = eps,
    forward_backward_ESS_threshold = forward_backward_ESS_threshold,
    method = method, n_max = n_max, n_threads = n_threads, smoother = smoother,
    Q_tilde = Q_tilde, est_a_0 = est_a_0, N_smooth_final = N_smooth_final,
    nu = nu)

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
    N_smooth_final <= N_smooth)
  stopifnot(
    typeof(nu) %in% c("double", "integer"), length(nu) == 1L, nu >= 0L,
    as.integer(nu) == nu)

  return(control)
}


.get_PF_static_args <- function(
  formula, data, by, max_T = NULL, id, trace, model, order, fixed = NULL,
  random = NULL){
  # get design matrix and risk set
  tmp <- get_design_matrix_and_risk_obj(
    formula = formula, data = data, by = by,
    max_T = if(is.null(max_T)) NULL else max_T, verbose = trace > 0,
    is_for_discrete_model = model == "logit", id = id, fixed = fixed,
    random = random)

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

  R_top <- t(R)
  out <- pf_fixed_effect_get_QR(
    clouds = clouds, risk_obj = risk_obj, ran_vars = X,
    fixed_terms = fixed_terms, R_top = R_top, tstart = tstart,
    tstop = tstop, fixed_parems = fixed_parems, family = family_arg,
    max_threads = nthreads, debug = debug)

  f_stack <- do.call(c, lapply(out, "[[", "f"))
  R_stack <- do.call(rbind, lapply(out, .get_R))

  qr. <- qr(R_stack, LAPACK = TRUE)
  f <- qr.qty(qr., f_stack)[1:nrow(fixed_terms)]

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
#' @export
get_Q_0 <- function(Qmat, Fmat){
  eg  <- eigen(Fmat)
  las <- eg$values
  if(any(abs(las) >= 1))
    stop("Divergent series")
  U   <- eg$vectors
  U_t <- t(U)
  T.  <- crossprod(U, Qmat %*% U)
  Z   <- T. / (1 - tcrossprod(las))
  out <- solve(U_t, t(solve(U_t, t(Z))))
  if(is.complex(out)){
    if(all(abs(Im(out)) < .Machine$double.eps^(3/4)))
      return(Re(out))

    stop("Q_0 has imaginary part")
  }

  out
}
