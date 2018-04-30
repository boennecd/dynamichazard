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
#' @param trace argument to get progress information. Zero will yield no info and larger integer values will yield incrementally more information.
#' @param model either \code{'logit'} for binary outcomes or \code{'exponential'} for piecewise constant exponential distributed arrival times.
#' @param seed seed to set at the start of every EM iteration.
#' @param ... optional way to pass arguments to \code{control}.
#'
#' @details
#' See \code{vignette("Particle_filtering", "dynamichazard")} for details.
#'
#' @return
#' An object of class \code{PF_EM}.
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
#'    n_threads = parallel::detectCores()),
#'  trace = 1)
#'
#'# Plot state vector estimates
#'plot(pf_fit, cov_index = 1)
#'plot(pf_fit, cov_index = 2)
#'plot(pf_fit, cov_index = 3)
#'
#'# Plot log-likelihood
#'plot(pf_fit$log_likes)}
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
#' @export
PF_EM <- function(
  formula, data, model = "logit", by, max_T, id, a_0, Q_0, Q, order = 1,
  control = PF_control(...), trace = 0, seed = .Random.seed, ...){
  #####
  # checks
  if(order != 1)
    stop(sQuote('order'), " not equal to 1 is not supported")

  if(!model %in% c("logit", "exponential"))
    stop(sQuote('model'), " is not supported")

  if(missing(id)){
    if(trace > 0)
      warning("You did not parse and Id argument")
    id = 1:nrow(data)
  }

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

  model_args <- .get_state_eq_matrices_PF(
    order = order, n_params = nrow(static_args$X),
    n_fixed = static_args$n_fixed, Q_0 = if(missing(Q_0)) NULL else Q_0,
    Q = if(missing(Q)) NULL else Q, a_0 = a_0)

  out <- do.call(.PF_EM, c(static_args, model_args, control, list(
    trace = trace, seed = seed, fixed_parems = fixed_parems)))

  out$call <- match.call()
  out
}

#' @importFrom graphics plot
.PF_EM <- function(
  n_fixed_terms_in_state_vec, X, fixed_terms, tstart, tstop, Q_0, Q, a_0, F.,
  L, R, m, risk_obj, n_max, n_threads, N_fw_n_bw, N_smooth, N_first, eps,
  forward_backward_ESS_threshold = NULL, debug = 0, trace,
  method = "AUX_normal_approx_w_particles", seed = NULL, smoother, model,
  fixed_parems){
  cl <- match.call()
  n_vars <- nrow(X)
  fit_call <- cl
  fit_call[[1]] <- as.name("PF_smooth")
  fit_call[["Q_tilde"]] <- bquote(diag(0, .(n_vars)))
  fit_call[["F"]] <- fit_call[["F."]]
  fit_call[c("eps", "seed", "F.", "trace")] <- NULL

  if(is.null(seed))
    seed <- .Random.seed

  log_likes <- rep(NA_real_, n_max)
  log_like <- log_like_max <- -Inf
  for(i in 1:n_max){
    if(trace > 0){
      if(i != 1)
        cat("\n\n\n")
      cat("#######################\nStarting EM iteration", i, "\n")

      cat("a_0 is:\n")
      print(eval(fit_call$a_0, environment()))

      if(length(fixed_parems) > 0){
        cat("Fixed parameters are:\n")
        print(fixed_parems)
      }

      cat("chol(Q) is:\n")
      print(chol(eval(fit_call$Q, environment())))
    }
    log_like_old <- log_like
    log_like_max <- max(log_like, log_like_max)

    #####
    # find clouds
    set.seed(seed)
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
    a_0_old <- eval(fit_call$a_0, environment())
    Q_old <- eval(fit_call$Q, environment())

    sum_stats <- compute_summary_stats(
      clouds, n_threads, a_0 = a_0, Q = Q, Q_0 = Q_0)
    a_0 <- drop(sum_stats[[1]]$E_xs)
    Q <- Reduce("+", lapply(sum_stats, "[[", "E_x_less_x_less_one_outers"))
    Q <- Q / length(sum_stats)

    fit_call$a_0 <- a_0
    fit_call$Q <- Q

    #####
    # Update fixed effects
    if(length(fixed_parems) > 0){
      if(trace > 0)
        cat("Updating fixed effects...\n")

      fit_call$fixed_parems <- fixed_parems <- .PF_update_fixed(
        clouds = clouds$smoothed_clouds, risk_obj = risk_obj, model = model,
        L = L, X = X, fixed_terms = fixed_terms, fixed_parems = fixed_parems,
        nthreads = n_threads)
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

    if(trace > 0)
      cat("The relative norm of the change in a_0 and Q are",
          a_0_relative_norm, "and", Q_relative_norm, "at iteration", i)

    if(has_converged <-
       Q_relative_norm < eps &&
       a_0_relative_norm < eps)
      break
  }

  if(!has_converged)
    warning("Method did not converge.")

  if(!exists("effective_sample_size", envir = environment()))
    effective_sample_size <- PF_effective_sample_size(clouds)

  return(structure(list(
    call = cl, clouds = clouds, a_0 = a_0, fixed_effects = fixed_parems, Q = Q,
    F = fit_call$F., L = L, R = R, summary_stats = sum_stats,
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
  method = "AUX_normal_approx_w_particles", n_max = 25,
  n_threads = getOption("ddhazard_max_threads"), smoother = "Fearnhead_O_N"){
  control <- list(
    N_fw_n_bw = N_fw_n_bw, N_smooth = N_smooth, N_first = N_first, eps = eps,
    forward_backward_ESS_threshold = forward_backward_ESS_threshold,
    method = method, n_max = n_max, n_threads = n_threads, smoother = smoother)

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

.get_state_eq_matrices_PF <- function(order, n_params, n_fixed, Q_0, Q, a_0){
  call. <- match.call()
  call.[["est_fixed_in_E"]] <- FALSE
  call.[[1]] <- quote(get_state_eq_matrices)
  out <- eval(call., envir = parent.frame())
  out$indicies_fix <- NULL

  out
}


.PF_update_fixed <- function(
  clouds, risk_obj,L, X, fixed_terms, fixed_parems, model, nthreads){
  # TODO: move this check
  if(!model %in% "logit")
    stop(sQuote(model), " is not implemented with fixed effects")

  family_arg <- switch (model, logit = "binomial")

  n_ps <- sapply(lapply(clouds, "[[", "states"), ncol)
  max_bit <- 400e6 # max byte to use on intermediate design matrix

  out <- NULL
  for(i in 1:length(clouds)){
    cl <- clouds[[i]]
    risk_set <- risk_obj$risk_sets[[i]]

    n_i <- length(risk_set)
    X_i <- X[, risk_set, drop = FALSE]
    X_i <- t(X_i) # note the transpose. Improves speed later
    fixed_terms_i <- fixed_terms[, risk_set, drop = FALSE]
    y_i <- risk_obj$is_event_in[risk_set] == (i - 1)

    n_ps_i <- n_ps[i]
    good <- which(drop(cl$weights) >= 1e-7)
    n_good <- length(good)

    # make `chunk` vectors and matrices to exploit parallel implementation
    chunk_size <- as.integer(max_bit /(8 * nrow(fixed_terms_i)))
    n_per_chunk <- min(as.integer((chunk_size + n_i) / n_i), length(good))
    fixed_terms_i_chunk <- matrix(
      # use that matrices are stored in column-major order
      rep(fixed_terms_i, n_per_chunk), nrow(fixed_terms_i),
      dimnames = list(dimnames(fixed_terms_i)[[1]], NULL))
    y_i_chunk <- rep(y_i, n_per_chunk)

    # Find weights and offsets and make QR decompositions
    particle_coefs <- L %*% cl$states
    i_start <- 1L
    while(i_start < n_good){
      i_end <- min(i_start + n_per_chunk - 1L, n_good)
      this_chunk <- good[i_start:i_end]

      offset. <- drop(X_i %*% particle_coefs[, this_chunk, drop = FALSE])
      weights. <- drop(rep_vec(cl$weights[this_chunk, ] * n_ps_i, n_i))

      if(length(this_chunk) == n_per_chunk){
        X_arg <- fixed_terms_i_chunk
        y_arg <- y_i_chunk

      } else {
        chunk_idx <- 1:(length(this_chunk) * n_i)
        X_arg <- fixed_terms_i_chunk[, chunk_idx, drop = FALSE]
        y_arg <- y_i_chunk[chunk_idx]

      }

      out <- c(out, list(
        parallelglm_QR_get_R_n_f(
          X = X_arg, Ys = y_arg, family = family_arg, beta0 = fixed_parems,
          weights = weights., offsets = offset., nthreads = nthreads)))

      i_start <- i_end + 1L
    }

    f_stack <- do.call(c, lapply(out, "[[", "f"))
    R_stack <- lapply(out, .get_R)
    R_stack <- do.call(rbind, R_stack)

    qr. <- qr(R_stack, LAPACK = TRUE)
    f <- qr.qty(qr., f_stack)[1:nrow(fixed_terms_i)]

    out <- list(list(
      R = qr.R(qr.), f = f,
      pivot = qr.$pivot - 1)) # less one to have zero index as cpp code
  }

  R <- .get_R(out[[1]])
  drop(solve(t(R) %*% R, t(R) %*% out[[1]]$f))
}

.get_R <- function(o){
  piv <- drop(o$pivot) + 1
  piv[piv] <- 1:length(piv)
  o$R[, piv, drop = FALSE]
}
