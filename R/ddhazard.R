if(getRversion() >= "2.15.1")
  utils::globalVariables(c("L", "R", "m", "indicies_fix"))

#' @title Fitting Dynamic Hazard Models
#' @description  Function to fit dynamic hazard models using state space models.
#' @param formula \code{\link[survival]{coxph}} like formula with \code{\link[survival]{Surv}(tstart, tstop, event)} on the left hand site of \code{~}.
#' @param data \code{data.frame} or environment containing the outcome and covariates.
#' @param model \code{"logit"}, \code{"cloglog"}, or \code{"exponential"} for
#' respectively the logistic link function with discrete outcomes,
#' the inverse cloglog link function with discrete outcomes,
#' or for the continuous time model with piecewise constant exponentially
#' distributed arrival times.
#' @param by interval length of the bins in which parameters are fixed.
#' @param max_T end of the last interval interval.
#' @param id vector of ids for each row of the in the design matrix.
#' @param a_0 vector \eqn{a_0} for the initial coefficient vector for the first iteration (optional). Default is estimates from static model (see \code{\link{static_glm}}).
#' @param Q_0 covariance matrix for the prior distribution.
#' @param Q initial covariance matrix for the state equation.
#' @param order order of the random walk.
#' @param weights weights to use if e.g. a skewed sample is used.
#' @param control list of control variables from
#' \code{\link{ddhazard_control}}.
#' @param verbose \code{TRUE} if you want status messages during execution.
#'
#' @details
#' This function can be used to estimate survival models where the regression parameters follows a given order random walk. The order is specified by the \code{order} argument. 1. and 2. order random walks is implemented. The regression parameters are updated at time \code{by}, 2\code{by}, ..., \code{max_T}. See the \code{vignette("ddhazard", "dynamichazard")} for details.
#'
#' All filter methods needs a state covariance matrix \code{Q_0} and state vector \code{a_0}. An estimate from a time-invariant model is used for \code{a_0} if it is not supplied (the same model you would get from \code{\link{static_glm}}). A diagonal matrix with large entries is recommended for \code{Q_0}. What is large dependents on the data set and \code{model}. Further, a covariance matrix for the first iteration \code{Q} is needed. The \code{Q} and \code{a_0} are estimated with an EM-algorithm.
#'
#' The model is specified through the \code{model} argument. The discrete
#' outcome models are where outcomes are binned into the intervals.
#' Be aware that there can be "loss" of information due to binning if outcomes
#' are not discrete to start with. It is key for these models that the
#' \code{id} argument is provided if individuals in the data set have
#' time-varying covariates. The the exponential model use a piecewise constant
#' exponential distribution
#' for the arrival times where there is no "loss" information due to binning.
#' Though, one of the assumptions of the model is not satisfied if outcomes
#' are only observed in discrete time intervals.
#'
#' It is recommended to see the Shiny app demo for this function by calling \code{\link{ddhazard_app}()}.
#'
#' @return
#' A list with class \code{ddhazard}. The list contains
#' \item{formula}{the passed formula.}
#' \item{call}{the matched call.}
#' \item{state_vecs}{2D matrix with the estimated state vectors (regression parameters) in each bin.}
#' \item{state_vars}{3D array with smoothed variance estimates for each state vector.}
#' \item{lag_one_cov}{3D array with lagged correlation matrix for each for each change in the state vector. Only present when the model is logit and the method is EKF.}
#' \item{n_risk}{the number of observations in each interval.}
#' \item{times}{the interval borders.}
#' \item{risk_set}{the object from \code{\link{get_risk_obj}} if saved.}
#' \item{data}{the \code{data} argument if saved.}
#' \item{weights}{\code{weights} used in estimation if saved.}
#' \item{id}{ids used to match rows in \code{data} to individuals.}
#' \item{order}{order of the random walk.}
#' \item{F_}{matrix which map from one state vector to the next.}
#' \item{method}{method used in the E-step.}
#' \item{est_Q_0}{\code{TRUE} if \code{Q_0} was estimated in the EM-algorithm.}
#' \item{family}{Rcpp \code{\link{Module}} with C++ functions used for estimation given the \code{model} argument.}
#' \item{discrete_hazard_func}{the hazard function corresponding to the \code{model} argument.}
#' \item{terms}{the \code{\link{terms}} object used.}
#' \item{has_fixed_intercept}{\code{TRUE} if the model has a time-invariant intercept.}
#' \item{xlev}{a record of the levels of the factors used in fitting.}
#'
#' @seealso
#' \code{\link[=plot.ddhazard]{plot}}, \code{\link[=residuals.ddhazard]{residuals}}, \code{\link[=predict.ddhazard]{predict}}, \code{\link{static_glm}}, \code{\link{ddhazard_app}}, \code{\link{ddhazard_boot}}
#'
#' @references
#' Fahrmeir, Ludwig. \emph{Dynamic modelling and penalized likelihood estimation for discrete time survival data}. Biometrika 81.2 (1994): 317-330.
#'
#' Durbin, James, and Siem Jan Koopman. \emph{Time series analysis by state space methods}. No. 38. Oxford University Press, 2012.
#'
#' @examples
#'# example with first order model
#'library(dynamichazard)
#'fit <- ddhazard(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
#'  Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
#'  control = ddhazard_control(method = "GMA"))
#'plot(fit)
#'
#'# example with second order model
#'fit <- ddhazard(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
#'  Q_0 = diag(1, 4), Q = diag(1e-4, 2), by = 50,
#'  control = ddhazard_control(method = "GMA"),
#'  order = 2)
#'plot(fit)
#'
#' @export
ddhazard = function(
  formula, data, model = "logit", by, max_T, id, a_0, Q_0, Q = Q_0, order = 1,
  weights, control = ddhazard_control(), verbose = F){
  #####
  # checks
  if (model %in% c("exp_bin", "exp_clip_time", "exp_clip_time_w_jump")){
    message(sQuote(model), " is not used after version 0.5.0.",
            " Use ", sQuote("exponential"), " instead.")
    model <- "exponential"

  }

  if(missing(id)){
    if(verbose)
      warning("You did not parse and Id argument. You should do so if you",
              " have time-varying covariates")
    id = 1:nrow(data)
  }

  if(missing(weights))
    weights <- rep(1, nrow(data))

  if(length(weights) != nrow(data)){
    stop("Length of weights does not match number of rows in data")

  } else if(!missing(id) && length(weights) != length(id)){
    stop("Length of weights does not match number of ids")

  } else{
    data <- data[weights > 0, ]
    id <- id[weights > 0]
    weights <- weights[weights > 0]

  }

  if(any(weights != 1))
    message("lag_one_cov will not be correct when some weights are not 1")

  if(model %in% c("logit", "cloglog")){
    is_for_discrete_model <- TRUE

  } else if (model %in% exp_model_names){
    is_for_discrete_model <- FALSE

  } else
    stop("Model ", sQuote(model), " is not implemented")

  #####
  # find risk set and design matrix
  tmp <- get_design_matrix_and_risk_obj(
    formula = formula, data = data, by = by,
    max_T = if(missing(max_T)) NULL else max_T, verbose = verbose,
    is_for_discrete_model = is_for_discrete_model, id = id)

  X_Y <- tmp$X_Y
  n_params <- tmp$n_params
  n_fixed <- tmp$n_fixed
  risk_set <- tmp$risk_set
  rm(tmp)

  #####
  # set control arguments
  ctrl <- substitute(control)
  if(is.call(ctrl) && identical(ctrl[[1]], quote(list))){
    # TODO: remove this warning sometime post version 0.6.0
    warning(sQuote("ddhazard_control"), " instead of ", sQuote("list"),
            " to the ", sQuote("control"), " argument")
    ctrl[[1L]] <- quote(ddhazard_control)
    control <- eval(ctrl)

  } else if(is.list(control)){
    if(any(!names(control) %in% c(
      names(formals(ddhazard_control)), names(ddhazard_control_xtra_get()))))
      stop(sQuote("control"), " contains unused elements")

    # does control have all the needed elements?
    if(!setequal(names(control), names(formals(ddhazard_control)))){
      cl <- as.call(c(list(quote(ddhazard_control)), control))
      control <- eval(cl)

    }

  } else
    stop("Invalid ", sQuote("control"), " argument")

  est_fixed_in_E <- control$fixed_terms_method == "E_step" && n_fixed > 0

  #####
  # Find starting values at time zero
  tmp <- get_start_values(
    data = data, formula = formula, max_T = max_T,
    # TODO: avoid transpose here by transpoing earlier
    X = t(X_Y$X), fixed_terms = t(X_Y$fixed_terms),
    risk_set = risk_set,
    verbose = verbose, n_threads = control$n_threads, model = model,
    a_0 = if(missing(a_0)) NULL else a_0,
    order = order, type = "RW",
    fixed_parems_start = control$fixed_parems_start)

  a_0 <- tmp$a_0
  control$fixed_parems_start <- fixed_parems_start <- tmp$fixed_parems_start

  #####
  # Find matrices for state equation
  tmp <- get_state_eq_matrices(
    order = order, n_params = n_params, n_fixed = n_fixed,
    est_fixed_in_E = est_fixed_in_E,
    Q_0 = if(missing(Q_0)) NULL else Q_0,
    Q = if(missing(Q)) NULL else Q,
    a_0 = a_0, control)
  list2env(tmp, environment())

  .check_filter_input(
    Q = Q, Q_0 = Q_0, F. = F., R = R, a_0 = a_0, L = L,
    fixed_parems = fixed_parems_start, est_fixed_in_E = est_fixed_in_E,
    X = X_Y$X, fixed_terms = X_Y$fixed_terms, order = order,
    has_transposed_design = FALSE)

  if(verbose)
    report_pre_liminary_stats_before_EM(risk_set = risk_set, Y = X_Y$Y)

  if(control$permu){
    # Permuting is useful for e.g. the SMA
    eval(get_permu_data_exp(X_Y[1:3], risk_set, weights))

  }

  if(!est_fixed_in_E){
    # we transpose for performance due to the column-major ordering. The logic
    # is that we primary look up indviduals (i.e. columns in the tranpose)
    X_Y$X <- t(X_Y$X)
    X_Y$fixed_terms <- t(X_Y$fixed_terms) # same

  } else {
    # We add the fixed covariates to the design matrix. Later, we recover the
    # fixed covariates
    X_Y$X <- t(cbind(X_Y$X, X_Y$fixed_terms))
    X_Y$fixed_terms <- matrix(nrow = 0, ncol = ncol(X_Y$X))

  }

  result <- NA
  k_vals <- seq_len(control$LR_max_try) - 1
  for(k in k_vals){
    tryCatch({
      result <- ddhazard_no_validation(
        a_0 = a_0, Q_0 = Q_0, F. = F., verbose = verbose, Q = Q,
        risk_set= risk_set, X_Y = X_Y, model = model, R = R, L = L,
        LR = control$LR * control$LR_decrease_fac^(k),
        n_fixed_terms_in_state_vec = ifelse(est_fixed_in_E, n_fixed, 0),
        weights = weights,
        control)
    }, error = function(e)
      if(length(k_vals) == 1 ||
          !grepl("^ddhazard_fit_cpp estimation error:", e$message))
        stop(e))

    LR <- control$LR * control$LR_decrease_fac^(k)
    if(did_fit <- !(length(result) == 1 && is.na(result)))
      break
  }

  # Check if fit succeeded
  if(!did_fit)
    stop("Failed to estimate model. The following learning rates was tried: ",
         paste0(control$LR * control$LR_decrease_fac^k_vals, collapse = ", "),
         ". Try decreasing the learning rate or change denom_term")

  if(control$method != "EKF")
    result$lag_one_cov = NULL # only works w/ EKF

  if(k != 0)
    message("Did not succed to fit the model wit a learning rate of ",
            control$LR, ". The learning rate was decrease by a factor ",
            control$LR_decrease_fac^k, " to yield a fit")

  if(est_fixed_in_E){
    # Recover the X and fixed term design matricies
    X_Y$fixed_terms <- X_Y$X[n_params + 1:n_fixed, , drop = F]
    if(n_params > 0)
      X_Y$X <- X_Y$X[1:n_params, , drop = F] else
        X_Y$X <- matrix(nrow = 0, ncol = ncol(X_Y$X))

    # We need to change the dimension of various arrays
    result$V_t_d_s <- result$V_t_d_s[-indicies_fix, -indicies_fix, , drop = F]

    result$fixed_effects <- result$a_t_d_s[1, indicies_fix]
    result$a_t_d_s <- result$a_t_d_s[, -indicies_fix, drop = F]

    result$B_s <- result$B_s[-indicies_fix, -indicies_fix, , drop = F]

    result$lag_one_cov <-
      result$lag_one_cov[-indicies_fix, -indicies_fix, , drop = F]

    F. <- F.[-indicies_fix, -indicies_fix, drop = F]
  }

  if(control$permu)
    eval(get_permu_data_rev_exp(list(), risk_set, weights))

  # Set names
  tmp_names = rep(rownames(X_Y$X), order)
  colnames(result$a_t_d_s) <- tmp_names
  dimnames(result$V_t_d_s) <- list(tmp_names, tmp_names, NULL)
  dimnames(result$Q) <- list(rownames(X_Y$X), rownames(X_Y$X))

  result$fixed_effects <- c(result$fixed_effects)
  names(result$fixed_effects) <- fixed_names <- gsub(
    "(^ddFixed\\()(.+)(\\)$)","\\2", rownames(X_Y$fixed_terms))

  if(est_fixed_in_E){
    .names <- if(order == 1)
      c(rownames(X_Y$X), fixed_names) else
        c(rownames(X_Y$X), fixed_names, rownames(X_Y$X))
    dimnames(result$Q_0) = list(.names, .names)

  } else
    dimnames(result$Q_0) <- list(tmp_names, tmp_names)

  if(model == "logit") {
    family <- Module("dd_logistic")

  } else if(model %in% exp_model_names) {
    family <- Module("dd_exponential")

  } else if(model == "cloglog") {
    family <- Module("dd_cloglog")

  } else {
    family <- NULL
    warning(sQuote("family"), " element not set on returned object")

  }

  discrete_hazard_func <- function(eta, at_risk_length){
    out <-  family$log_like(outcome = logical(length(eta)), eta = eta,
                            at_risk_length = at_risk_length)
    1 - exp(out)
  }

  structure(list(
    family = family,
    discrete_hazard_func = discrete_hazard_func,
    formula = formula,
    call = match.call(),
    state_vecs = result$a_t_d_s,
    state_vars = result$V_t_d_s,
    lag_one_cov = result$lag_one_cov,
    fixed_effects = result$fixed_effects,
    n_iter = result$n_iter,
    Q = result$Q,
    Q_0 = result$Q_0,
    n_risk = unlist(lapply(risk_set$risk_sets, length)),
    times = risk_set$event_times,
    risk_set = if(control$save_risk_set) risk_set else NULL,
    data = if(control$save_data) data else NULL,
    weights = if(control$save_data) weights else NULL,
    id = if(control$save_data) id else NULL,
    order = order, F_ = F.,
    method = control$method,
    model = model,
    est_Q_0 = control$est_Q_0,
    terms = X_Y$terms,
    has_fixed_intercept = X_Y$has_fixed_intercept,
    xlev = X_Y$xlev,
    control = control,
    LR = LR),
    "class" = "ddhazard")
}

ddhazard_no_validation <- function(
  a_0, Q_0, F., verbose, Q, risk_set, X_Y, model, LR, n_fixed_terms_in_state_vec,
  weights = weights, control, R, L){
  ddhazard_fit_cpp(
    a_0 = a_0, Q_0 = Q_0, F_ = F., verbose = verbose,
    Q = Q, n_max = control$n_max,
    R = R, L = L,
    risk_obj = risk_set, eps = control$eps,
    X = X_Y$X, fixed_terms = X_Y$fixed_terms,
    fixed_parems_start = control$fixed_parems_start,
    tstart = X_Y$Y[, 1], tstop = X_Y$Y[, 2],
    est_Q_0 = control$est_Q_0, method = control$method,
    model = model,
    kappa = control$kappa, alpha = control$alpha, beta = control$beta,
    NR_eps = control$NR_eps,
    LR = LR,
    eps_fixed_parems = control$eps_fixed_parems,
    max_it_fixed_params = control$max_it_fixed_params,
    fixed_effect_chunk_size = control$fixed_effect_chunk_size,
    debug = control$debug,
    n_threads = control$n_threads,
    denom_term = control$denom_term,
    n_fixed_terms_in_state_vec = n_fixed_terms_in_state_vec,
    weights = weights,
    use_pinv = control$use_pinv,
    criteria = control$criteria,
    posterior_version = control$posterior_version,
    GMA_max_rep = control$GMA_max_rep,
    GMA_NR_eps = control$GMA_NR_eps,
    EKF_batch_size = control$EKF_batch_size,
    est_a_0 = control$est_a_0)
}

get_state_eq_matrices <-  function(
  order, n_params, n_fixed, est_fixed_in_E = FALSE, Q_0, Q, a_0, control,
  type = "RW", F. = NULL){
  #####
  # Indices vector used later
  state_dim <- n_params * order + n_fixed * est_fixed_in_E
  indicies_cur <- 1:n_params
  indicies_lag <- if(order == 2)
    n_params + 1:n_params + ifelse(est_fixed_in_E, n_fixed, 0) else vector()
  indicies_fix <- if(est_fixed_in_E) 1:n_fixed + n_params else vector()

  #####
  # Checks
  func <- function(x, n_elem, default){
    if(!is.null(x))
      return(x)

    if(n_elem > 0)
      warning(
        sQuote(substitute(x)),
        " not supplied. It has been set to a diagonal matrix with diagonal entries equal to ",
        default)

    diag(default, n_elem)
  }

  Q_0 <- func(Q_0, n_params * order, 10)
  Q <- func(Q, n_params, 1)

  if(is.vector(Q) && length(Q) == 1)
    Q <- matrix(Q)

  if(is.vector(Q_0) && length(Q_0) == 1)
    Q_0 <- matrix(Q_0)

  if(!order %in% 1:2)
    stop("Method not implemented for order ", order)

  #####
  # Get F matrix
  if(is.null(F.)){
    if(type == "RW"){
      if(order == 1){
        F. <- diag(rep(1, n_params + n_fixed * est_fixed_in_E))

      } else if(order == 2){
        F. = matrix(NA_real_,
                    nrow = 2 * n_params + n_fixed * est_fixed_in_E,
                    ncol = 2 * n_params + n_fixed * est_fixed_in_E)
        F.[indicies_cur, indicies_cur] = diag(2, n_params)
        F.[indicies_lag, indicies_cur] = diag(1, n_params)
        F.[indicies_cur, indicies_lag] = diag(-1, n_params)
        F.[indicies_lag, indicies_lag] = 0

        if(length(indicies_fix) > 0){
          F.[indicies_fix, ] <- 0
          F.[, indicies_fix] <- 0
          if(length(indicies_fix) > 1)
            diag(F.[indicies_fix, indicies_fix]) <- 1 else
              F.[indicies_fix, indicies_fix] <- 1
        }
      }
    } else if(type == "VAR" && order == 1 && !est_fixed_in_E){
      F. <- diag(.1, n_params) # TODO: find better default
    } else
      stop("Default parameter not implemented for ", sQuote("F."))
  }

  #####
  # Setup for Q_0 and a_0
  if(est_fixed_in_E){
    # We need to add entries to the matrix and vector

    if(ncol(Q_0) != n_params * order + length(indicies_fix)){
      Q_0_new <- F.
      Q_0_new[, ] <- 0
      Q_0_new[-indicies_fix, -indicies_fix] <- Q_0
      if(length(indicies_fix) == 1){
        Q_0_new[indicies_fix, indicies_fix] <- control$Q_0_term_for_fixed_E_step
      } else
        diag(Q_0_new[indicies_fix, indicies_fix]) <-
          control$Q_0_term_for_fixed_E_step
      Q_0 <- Q_0_new
    }

    if(order == 1){
      a_0 <- c(a_0, control$fixed_parems_start)

    } else if(order == 2){
      first_half <- 1:(length(a_0)/2)
      a_0 <- c(a_0[first_half], control$fixed_parems_start, a_0[-first_half])

    } else
      stop("Order of ", order, " not supported")
  }

  #####
  # Setup for other matrices and vectors
  R <- matrix(0, state_dim, n_params)
  if(n_params > 0)
    diag(R)[1:n_params] <- 1
  ncol_L <- n_params +  n_fixed * est_fixed_in_E
  L <- matrix(0, ncol_L, state_dim)
  if(n_params + n_fixed * est_fixed_in_E > 0)
    diag(L)[1:ncol_L] <- 1

  return(list(
    Q = Q, Q_0 = Q_0, F. = F., R = R, L = L, a_0 = a_0,
    indicies_fix = indicies_fix))
}

.check_filter_input <- function(
  Q, Q_0, F., R, a_0, L = NULL, fixed_parems, est_fixed_in_E,
  X, fixed_terms, order, has_transposed_design = TRUE, Q_tilde = NULL,
  G = NULL, J = NULL, K = NULL, theta = NULL, psi = NULL, phi = NULL){
  lp_dim  <- if(has_transposed_design) nrow(          X) else ncol(          X)
  fix_dim <- if(has_transposed_design) nrow(fixed_terms) else ncol(fixed_terms)
  state_dim <- lp_dim * order + fix_dim * est_fixed_in_E
  rng_dim   <- lp_dim

  .check_full_rank_square(Q      , rng_dim  , TRUE , is_null_ok = TRUE)
  .check_full_rank_square(Q_0    , state_dim, TRUE , is_null_ok = TRUE)
  .check_full_rank_square(F.     , state_dim, FALSE, is_null_ok = TRUE)
  .check_full_rank_square(Q_tilde, rng_dim  , TRUE , is_null_ok = TRUE,
                          is_zero_ok = TRUE)

  .check_selection_matrix(R, state_dim, rng_dim)
  if(!is.null(L))
    .check_selection_matrix(L, lp_dim + fix_dim * est_fixed_in_E, state_dim)

  if(!length(a_0) == state_dim)
    stop("Invalid ", sQuote("a_0"))
  if(!length(fixed_parems) == fix_dim)
    stop("Invalid ", sQuote("fixed_terms"))

  if(!is.null(G) || !is.null(theta)){
    stopifnot(nrow(G) == state_dim^2, ncol(G) == length(theta),
              length(theta) <= nrow(G), qr(G)$rank == ncol(G))
    F. <- .get_F(G, theta)
    if(qr(F.)$rank < ncol(F.))
      stop(sQuote("F"), " is singular with the starting values")

  }
  if(!is.null(J) || !is.null(psi) || !is.null(K) || !is.null(phi)){
    stopifnot(nrow(J) == rng_dim, ncol(J) == length(psi),
              length(psi) <= nrow(J), qr(J)$rank == ncol(J))
    stopifnot(nrow(K) == rng_dim * (rng_dim - 1) / 2, ncol(K) == length(phi),
              length(phi) <= nrow(K), qr(K)$rank == ncol(K))

    Q <- .get_Q(J, K, psi, phi)$Q
    eg <- eigen(Q)

    if(any(abs(eg$values) < 1e-10))
      stop(sQuote("Q"), " is not positive definite with the starting values")
  }


  invisible(TRUE)
}

.check_full_rank_square <- function(X, expected_dim, pos_def,
                                    is_null_ok = FALSE, is_zero_ok = FALSE){
  qu <- substitute({
    if(!is.null(X) || !is_null_ok){
      if(ncol(X) != n || nrow(X) != n)
        stop("Invalid dimensions of ", sQuote(Xstr), ". Should be (", n,
             ", ", n, ") but is ",
             paste0("(", paste0(dim(X), collapse = ", "), ")"))
      if(n > 0){
        if(pos_def){
          eg <- eigen(X)
          tol <- if(is_zero_ok) 0. else
            .Machine$double.eps * NCOL(X) * max(abs(eg$values))
          if(!all(eg$values >= tol))
            stop(sQuote(Xstr), " is not positive definite")
        } else
          if(qr(X)$rank < n)
            stop(sQuote(Xstr), " does not have full rank")
      }
    }
  }, list(X = substitute(X), Xstr = deparse(substitute(X)), n = expected_dim,
          pos_def = pos_def, is_null_ok = is_null_ok, is_zero_ok = is_zero_ok))

  eval(qu, envir = parent.frame())
}

.check_selection_matrix <- function(X, n, m){
  qu <- substitute({
    if(nrow(X) != n || ncol(X) != m)
      stop("Invalid dimensions of ", sQuote(Xstr),  ". Should be (", m,
           ", ", n, ") but is ",
           paste0("(", paste0(dim(X), collapse = ", "), ")"))
    if(!all(X %in% c(0, 1)))
      stop("All entries of ", sQuote(Xstr), " are not zero or one")
    if(n * m > 0 && !qr(X)$rank == min(n, m))
      stop(sQuote(Xstr), " rank is less than ", min(n, m))
  }, list(X = substitute(X), Xstr = deparse(substitute(X)), n = n, m = m))

  eval(qu, envir = parent.frame())
}

# TODO: remove other exp_ names at some future point after 0.5.0 changes
exp_model_names <- c(
  "exponential", "exp_bin", "exp_clip_time", "exp_clip_time_w_jump")

#' @importFrom utils capture.output
get_start_values <- function(
  data, formula, max_T, X, risk_obj, fixed_terms, risk_set, verbose = FALSE,
  n_threads, model, order, a_0 = NULL, fixed_parems_start = NULL,
  fixed = NULL, random = NULL, type){
  n_params = nrow(X)
  n_fixed = nrow(fixed_terms)

  miss_a_0 <- is.null(a_0) && n_params > 0
  miss_fix <- is.null(fixed_parems_start) && n_fixed > 0

  finish <- function(a = a_0, f = fixed_parems_start,
                     adj_a_0 = TRUE){
    if(adj_a_0)
      a <- rep(a, order)

    if(verbose){
      if(length(a) > 0)
        message("Starting values for time-varying coeffecients are:\n",
                paste0(capture.output(a), collapse = "\n"))
      if(length(f > 0))
        message("Starting values for fixed coeffecients are:\n",
                paste0(capture.output(f), collapse = "\n"))
    }

    list(a_0 = a, fixed_parems_start = f)
  }

  if(!miss_a_0 && !miss_fix){
    if(is.null(fixed_parems_start))
      fixed_parems_start <- numeric()
    if(is.null(a_0))
      a_0 <- numeric()

    return(finish(adj_a_0 = FALSE))

  }

  if(n_params == 0){
    # Model is fitted for testing
    warning("The model can be estimated more effeciently by using get_survival_case_Weights_and_data and static_glm when there is no time-varying parameters")
    a_0 = vector()

    if(is.null(fixed_parems_start))
      fixed_parems_start <- rep(0, n_fixed) else
        fixed_parems_start <- fixed_parems_start

    return(finish(adj_a_0 = TRUE))

  }

  #####
  # setup call with to use with different design matrices and families
  cl <- bquote(static_glm(
    formula = .(if(missing(formula)) quote(fixed) else quote(formula)),
    data = data, max_T = max_T, risk_obj = risk_set,
    epsilon = sqrt(.Machine$double.eps) * 1e1,
    speedglm = FALSE, only_coef = TRUE, mf = mt,
    method_use = "parallelglm_quick", n_threads = n_threads))

  if(model == "logit"){
    cl$family <- "binomial"

  } else if(model %in% exp_model_names){
    cl$family <- "exponential"

  } else if(model == "cloglog") {
    cl$family <- "cloglog"

  } else
    stop("Method not implemented to find initial values for ", sQuote(model),
         ". Please, provide intial values for a_0")

  if(type == "VAR"){
    if(miss_fix){
      mt <- t(fixed_terms)
      fixed_parems_start <- eval(cl)

    } else if(nrow(fixed_terms) == 0L)
      fixed_parems_start <- numeric()

    if(miss_a_0)
      a_0 <- numeric(nrow(X)) else if (nrow(X) == 0L)
        a_0 <- numeric()

  } else if (type == "RW") {
    mt <- cbind(t(X), t(fixed_terms))
    out <- eval(cl)
    if(miss_fix)
      fixed_parems_start <- out[-(1:nrow(X))] else if(nrow(fixed_terms) == 0L)
        fixed_parems_start <- numeric()

    if(miss_a_0)
      a_0 <- out[1:nrow(X)] else if(nrow(X) == 0L)
        a_0 <- numeric()

  } else
    stop("type ", sQuote(type), " is not supported")

  if(!miss_a_0 && miss_fix)
    return(finish(adj_a_0 = TRUE))

  # give an additional message
  message("a_0 not supplied. IWLS estimates of static glm model is used",
          " for random walk models. Otherwise the values are zero")
  return(finish(adj_a_0 = TRUE))
}

report_pre_liminary_stats_before_EM <- function(risk_set, Y){
  tmp_tbl = matrix(NA_real_, nrow = risk_set$d, ncol = 2)
  colnames(tmp_tbl) = c("Risk size", "Num events")

  # Find the number of events in each bin
  n_events <- xtabs(~ risk_set$is_event_in)
  n_events <- n_events[names(n_events) != "-1"]
  names(n_events) <- as.integer(names(n_events)) + 1

  # Some bins may have no events. We need to add these
  if(any(event_group_missing <- !seq_len(risk_set$d) %in% names(n_events))){
    n_missing <- sum(event_group_missing)
    n_events <- c(rep(0, n_missing), n_events)
    names(n_events)[seq_len(n_missing)] <- which(event_group_missing)
  }
  n_events <- n_events[order(as.integer(names(n_events)))]

  # Insert into table to print and print
  tmp_tbl[, "Num events"] <- n_events
  tmp_tbl[, "Risk size"] <- unlist(lapply(risk_set$risk_sets, length))

  message("Total number of included events are ", sum(tmp_tbl[, 2]), " of ",   sum(Y[, 3]))
  message("Size of risk set and number of events in each risk set are ([Risk set size]:[# events]):")

  tmp_tbl[, 1] <- sprintf("%8s", tmp_tbl[, 1])
  tmp_tbl[, 2] <- sprintf("%-8s", tmp_tbl[, 2])
  tmp_message <- sprintf("%21s", apply(tmp_tbl, 1, paste, collapse = " : "))
  msg_final <- tmp_message[1]
  for(i in seq_along(tmp_message)[-1])
    msg_final <- paste0(msg_final, if((i - 1) %% 4 > 0) " " else "\n", tmp_message[i])
  message(msg_final)

  message("Running EM")

  invisible()
}

get_design_matrix_and_risk_obj <- function(
  formula, data, by, max_T = NULL, id, verbose = FALSE, is_for_discrete_model,
  fixed = NULL, random = NULL){
  X_Y = get_design_matrix(formula, data, fixed = fixed, random = random)

  if(verbose)
    message("Finding Risk set")

  risk_set <-
    get_risk_obj(
      Y = X_Y$Y, by = by,
      max_T = ifelse(
        is.null(max_T),
        min(max(X_Y$Y[X_Y$Y[, 3] == 1, 2]), max(X_Y$Y[X_Y$Y[, 3] == 0, 2])),
        max_T),
      id = id, is_for_discrete_model = is_for_discrete_model)

  list(X_Y = X_Y, risk_set = risk_set,
       n_params = ncol(X_Y$X), n_fixed = ncol(X_Y$fixed_terms))
}

#' @title Auxiliary for Controlling Dynamic Hazard Models
#'
#' @description
#' Auxiliary for additional settings with \code{\link{ddhazard}}.
#'
#' @param method set to the method to use in the E-step. Either \code{"EKF"}
#' for the Extended Kalman Filter, \code{"UKF"} for the Unscented Kalman
#' Filter, \code{"SMA"} for the sequential posterior mode approximation method
#' or \code{"GMA"} for the global mode approximation method. \code{"EKF"} is
#' the default.
#' @param LR learning rate.
#' @param NR_eps tolerance for the Extended Kalman filter. Default is
#' \code{NULL} which means that no extra iteration is made in the correction
#' step.
#' @param alpha hyper parameter \eqn{\alpha} in the unscented Kalman Filter.
#' @param beta hyper parameter \eqn{\beta} in the unscented Kalman Filter.
#' @param kappa hyper parameter \eqn{\kappa} in the unscented Kalman Filter.
#' @param n_max maximum number of iteration in the EM-algorithm.
#' @param eps tolerance parameter for the EM-algorithm
#' @param est_Q_0 \code{TRUE} if you want the EM-algorithm to estimate
#' \code{Q_0}. Default is \code{FALSE}.
#' @param save_risk_set \code{TRUE} if you want to save the list from
#' \code{\link{get_risk_obj}} used to estimate the model. It may be needed for
#' later calls to e.g., \code{residuals}, \code{plot} and \code{logLike}.
#' @param save_data \code{TRUE} if you want to keep the \code{data} argument.
#' It may be needed for later calls to e.g., \code{residuals}, \code{plot} and
#' \code{logLike}.
#' @param denom_term term added to denominators in either the EKF or UKF.
#' @param n_threads maximum number of threads to use.
#' @param fixed_parems_start starting value for fixed terms.
#' @param fixed_terms_method the method used to estimate the fixed effects.
#' Either \code{'M_step'} or \code{'E_step'} for estimation in the M-step or
#' E-step respectively.
#' @param Q_0_term_for_fixed_E_step the diagonal value of the initial
#' covariance matrix, \code{Q_0}, for the fixed effects if fixed effects are
#' estimated in the E-step.
#' @param eps_fixed_parems tolerance used in the M-step of the Fisher's scoring
#' algorithm for the fixed effects
#' @param permu \code{TRUE} if the risk sets should be permutated before
#' computation. This is \code{TRUE} by default for posterior mode approximation
#' method and \code{FALSE} for all other methods.
#' @param posterior_version the implementation version of the posterior
#' approximation method. Either \code{"woodbury"} or \code{"cholesky"}.
#' @param GMA_max_rep maximum number of iterations in the correction step if
#' \code{method = 'GMA'}.
#' @param GMA_NR_eps tolerance for the convergence criteria for the relative
#' change in the norm of the coefficients in the correction step if
#' \code{method = 'GMA'}.
#' @param est_a_0 \code{FALSE} if the starting value of the state model should
#' be fixed.
#' @param ... additional undocumented arguments.
#'
#' @return
#' A list with components named as the arguments.
#'
#' @seealso \code{\link{ddhazard}}
#'
#' @export
ddhazard_control <- function(
  kappa = NULL, alpha = 1, beta = 0,
  NR_eps = NULL, LR = 1, n_max = 10^2, eps = 1e-3,
  est_Q_0 = F, method = "EKF", save_risk_set = T,
  save_data = T, eps_fixed_parems = 1e-4,
  fixed_parems_start = NULL,
  n_threads = getOption("ddhazard_max_threads"),
  denom_term = 1e-5,
  fixed_terms_method = "E_step",
  Q_0_term_for_fixed_E_step = NULL,
  permu = if(!is.null(method)) method == "SMA" else F,
  posterior_version = "cholesky", GMA_max_rep = 25, GMA_NR_eps = 1e-4,
  est_a_0 = TRUE, ...) {
  dts <- list(...)
  xtra <- ddhazard_control_xtra_get()
  ma <- match(names(dts), names(xtra), NA_integer_)
  if(anyNA(ma))
    stop("Unused arguments passed to ", sQuote("ddhazard_control"))
  if(length(ma) > 0)
    xtra[ma] <- dts

  control <- c(list(
    kappa = kappa, alpha = alpha, beta = beta, NR_eps = NR_eps, LR = LR,
    n_max = n_max, eps = eps, est_Q_0 = est_Q_0, method = method,
    save_risk_set = save_risk_set, save_data = save_data,
    eps_fixed_parems = eps_fixed_parems,
    fixed_parems_start = fixed_parems_start, n_threads = n_threads,
    denom_term = denom_term, fixed_terms_method = fixed_terms_method,
    Q_0_term_for_fixed_E_step = Q_0_term_for_fixed_E_step, permu = permu,
    posterior_version = posterior_version, GMA_max_rep = GMA_max_rep,
    GMA_NR_eps = GMA_NR_eps, est_a_0 = est_a_0), xtra)
  if(!control$criteria %in% c("delta_coef", "delta_likeli"))
    stop("Convergence criteria ", control$criteria, " not implemented")

  if(is.null(control$Q_0_term_for_fixed_E_step))
    control$Q_0_term_for_fixed_E_step <- ifelse(
      # quite abritary values (especially the former - the latter is not too
      # important).
      control$method %in% c("UKF", "GMA"), 1, 1e5)

  if(!control$fixed_terms_method %in% c("M_step", "E_step"))
    stop("fixed_terms_method method ", sQuote(control$fixed_terms_method),
         " is not implemented")

  if(control$denom_term <= 0){
    stop("Method not implemented with penalty term ",
         sQuote("control$denom_term"), " equal to ", control$denom_term)

  }

  control
}

ddhazard_control_xtra_get <- function()
  list(
    max_it_fixed_params = 25, fixed_effect_chunk_size = 1e4, debug = F,
    LR_max_try = 10, LR_decrease_fac = 0.9, use_pinv = FALSE,
    criteria = "delta_coef", EKF_batch_size = 500L)
