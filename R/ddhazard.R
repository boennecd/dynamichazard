#' @title Function to fit dynamic discrete hazard models
#' @description  Function to fit dynamic discrete hazard models using state space models
#' @param formula \code{\link[survival]{coxph}} like formula with \code{\link[survival]{Surv}(tstart, tstop, event)} on the left hand site of \code{~}
#' @param data Data frame or environment containing the outcome and co-variates
#' @param model \code{"logit"}, \code{"exp_clip_time_w_jump"}, \code{"exp_clip_time"} or \code{"exp_bin"} for the discrete time function using the logistic link function in the first case or for the continuous time model with different estimation method in the three latter cases (see the ddhazard vignette for details of the methods)
#' @param by Interval length of the bins in which parameters are fixed
#' @param max_T End of the last interval. The last stop time with an event is selected if the parameter is omitted
#' @param id Vector of ids for each row of the in the design matrix
#' @param a_0 Vector \eqn{a_0} for the initial coefficient vector for the first iteration (optional). Default is estimates from static model (see \code{\link{static_glm}})
#' @param Q_0 Covariance matrix for the prior distribution
#' @param Q Initial covariance matrix for the state equation
#' @param order Order of the random walk
#' @param weights Weights to use if e.g. a skewed sample is used
#' @param control List of control variables (see details below)
#' @param verbose \code{TRUE} if you want status messages during execution
#'
#' @details
#' This function can be used to estimate a binary regression where the regression parameters follows a given order random walk. The order is specified by the \code{order} argument. 1. and 2. order random walks is implemented. The regression parameters are updated at time \code{by}, 2\code{by}, ..., \code{max_T}. See the vignette 'ddhazard' for more details
#'
#' All filter methods needs a state covariance matrix \code{Q_0} and state vector \code{a_0}. An estimate from a time-invariant model is provided for \code{a_0} if it is not supplied (the same model you would get from \code{\link{static_glm}} function). A diagonal matrix with large entries is recommended for \code{Q_0}. What is large dependents on the data set and \code{model}. Further, a variance matrix for the first iteration \code{Q} is needed. It is recommended to select diagonal matrix with low values for the latter. The \code{Q}, \code{a_0} and optionally \code{Q_0} is estimated with an EM-algorithm
#'
#' The model is specified through the \code{model} argument. See the \code{model} in the argument above for details. The logistic model is where outcomes are binned into the intervals. Be aware that there can be loss of information due to binning. It is key for the logit model that the \code{id} argument is provided if individuals in the data set have time varying co-variates. The the exponential models use an exponential model for the arrival times where there is no loss information due to binning
#'
#' It is recommended to see the Shiny app demo for this function by calling \code{\link{ddhazard_app}()}
#'
#' @section Control:
#' The \code{control} argument allows you to pass a \code{list} to select additional parameters. See the vignette 'ddhazard' for more information on hyper parameters. Unspecified elements of the list will yield default values
#' \describe{
#' \item{\code{method}}{Set to the method to use in the E-step. Either \code{"EKF"} for the Extended Kalman Filter, \code{"UKF"}for the Unscented Kalman Filter, \code{"SMA"} for the squential posterior mode approximation method or \code{"GMA"} for the global mode approximation method. \code{"EKF"} is the default}
#' \item{\code{LR}}{Learning rate for the Extended Kalman filter}
#' \item{\code{NR_eps}}{Tolerance for the Extended Kalman filter. Default is \code{NULL} which means that no extra iteration is made in the correction step}
#' \item{\code{alpha}}{Hyper parameter \eqn{\alpha} in the Unscented Kalman Filter}
#' \item{\code{beta}}{Hyper parameter \eqn{\beta} in the Unscented Kalman Filter }
#' \item{\code{kappa}}{Hyper parameter \eqn{\kappa} in the Unscented Kalman Filter}
#' \item{\code{n_max}}{Maximum number of iteration in the EM-algorithm}
#' \item{\code{eps}}{Tolerance parameter for the EM-algorithm}
#' \item{\code{est_Q_0}}{\code{TRUE} if you want the EM-algorithm to estimate \code{Q_0}. Default is \code{FALSE}}
#' \item{\code{save_risk_set}}{\code{TRUE} if you want to save the list from \code{\link{get_risk_obj}} used to estimate the model. It may be needed for later call to \code{residuals}, \code{plot} and \code{logLike}. Can be set to \code{FALSE} to save memory}
#' \item{\code{save_data}}{\code{TRUE} if you want to save the list \code{data} argument. It may be needed for later call to \code{residuals}, \code{plot} and \code{logLike}. Can be set to \code{FALSE} to save memory}
#' \item{\code{denom_term}}{Term added to denominators in either the EKF or UKF}
#' \item{\code{fixed_parems_start}}{Starting value for fixed terms}
#' \item{\code{fixed_terms_method}}{The method used to estimate the fixed effects. Either \code{'M_step'} or \code{'E_step'} for estimation in the M-step or E-step respectively}
#' \item{\code{Q_0_term_for_fixed_E_step}}{The diagonal value of the initial covariance matrix, \code{Q_0}, for the fixed effects if fixed effects are estimated in the E-step}
#' \item{\code{eps_fixed_parems}}{Tolerance used in the M-step of the Fisher's Scoring Algorithm for the fixed effects}
#' \item{\code{permu}}{\code{TRUE} if the risk sets should be permutated before computation. This is \code{TRUE} by default for posterior mode approximation method and \code{FALSE} for all other methods}
#' \item{\code{posterior_version}}{The implementation version of the posterior approximation method. Either \code{"woodbury"} or \code{"cholesky"}}
#' \item{\code{GMA_max_rep}}{Maximum number of iterations in the correction step if \code{method = 'GMA'}}
#' \item{\code{GMA_NR_eps}}{Tolerance for the convergence criteria for the relative change in the norm of the coefficients in the correction step if \code{method = 'GMA'}}
#'}
#'
#' @return
#' A list with class \code{fahrmeier_94}. The list contains:
#' \describe{
#' \item{\code{formula}}{The passed formula }
#' \item{\code{state_vecs}}{2D matrix with the estimated state vectors (regression parameters) in each bin }
#' \item{\code{state_vars}}{3D array with smoothed variance estimates for each state vector }
#' \item{\code{lag_one_cov}}{3D array with lagged correlation matrix for each for each change in the state vector. Only present when the model is logit and the method is EKF }
#' \item{\code{n_risk}}{The number of observations in each interval }
#' \item{\code{times}}{The interval borders }
#' \item{\code{risk_set}}{The object from \code{\link{get_risk_obj}} if saved }
#' \item{\code{data}}{The \code{data} argument if saved }
#' \item{\code{id}}{ids used to match rows in \code{data} to individuals }
#' \item{\code{order}}{Order of the random walk }
#' \item{\code{F_}}{Matrix with that map transition from one state vector to the next }
#' \item{\code{method}}{Method used in the E-step}
#' \item{\code{est_Q_0}}{\code{TRUE} if \code{Q_0} was estimated in the EM-algorithm }
#' \item{\code{hazard_func}}{Hazard function }
#' \item{\code{hazard_first_deriv}}{First derivative of the hazard function with respect to the linear predictor}
#'}
#'
#' @seealso
#' \code{\link[=plot.fahrmeier_94]{plot}}, \code{\link[=residuals.fahrmeier_94]{residuals}}, \code{\link[=predict.fahrmeier_94]{predict}}, \code{\link{static_glm}}, \code{\link{ddhazard_app}}, \code{\link{ddhazard_boot}}
#'
#' @references
#' Fahrmeir, Ludwig. \emph{Dynamic modelling and penalized likelihood estimation for discrete time survival data}. Biometrika 81.2 (1994): 317-330.
#'
#' Durbin, James, and Siem Jan Koopman. \emph{Time series analysis by state space methods}. No. 38. Oxford University Press, 2012.
#'
#' @export
ddhazard = function(formula, data,
                    model = "logit",
                    by, max_T, id,
                    a_0, Q_0, Q = Q_0,
                    order = 1, weights,
                    control = list(),
                    verbose = F){

  if(model == "exp_trunc_time"){
    message("Model 'exp_trunc_time' have been renamed to 'exp_clip_time'")
    model <- "exp_clip_time"
  } else if(model == "exp_trunc_time_w_jump"){
    message("Model 'exp_trunc_time_w_jump' have been renamed to 'exp_clip_time_w_jump'")
    model <- "exp_clip_time_w_jump"
  } else if(model == "exp_combined"){
    stop("'exp_combined' is not supported since version 0.3.0")
  }

  if(missing(id)){
    if(verbose)
      warning("You did not parse and Id argument")
    id = 1:nrow(data)
  }

  if(missing(weights)){
    weights <- rep(1, nrow(data))
  }

  if(length(weights) != nrow(data)){
    stop("Length of weights does not match number of rows in data")
  } else if(!missing(id) && length(weights) != length(id)){
    stop("Length of weights does not match number of ids")
  } else{
    data <- data[weights > 0, ]
    id <- id[weights > 0]
    weights <- weights[weights > 0]
  }

  if(any(weights != 1)){
    message("lag_one_cov will not be correct when some weights are not 1")
  }

  X_Y = get_design_matrix(formula, data)
  n_params = ncol(X_Y$X)

  if(model == "logit"){
    is_for_discrete_model <- TRUE

  } else if (model %in% exp_model_names){
    is_for_discrete_model <- FALSE

  } else
    stop("Model '", model, "' is not implemented")

  control_default <- list(kappa = NULL, alpha = 1, beta = 0,
                          NR_eps = NULL, LR = 1, n_max = 10^2, eps = .01,
                          est_Q_0 = F, method = "EKF", save_risk_set = T,
                          save_data = T, eps_fixed_parems = 1e-3,
                          max_it_fixed_params = 10, fixed_effect_chunk_size = 1e4,
                          debug = F, fixed_parems_start = NULL, LR_max_try = 10,
                          LR_decrease_fac = 0.9,
                          n_threads = getOption("ddhazard_max_threads"),
                          denom_term = .0001,
                          fixed_terms_method = "E_step",
                          Q_0_term_for_fixed_E_step = NULL,
                          use_pinv = T, criteria = "delta_coef",
                          permu = if(!is.null(control$method))
                            control$method == "SMA" else F,
                          posterior_version = "cholesky",
                          GMA_max_rep = 10,
                          GMA_NR_eps = 0.1,
                          EKF_batch_size = 2000L)

  if(any(is.na(control_match <- match(names(control), names(control_default)))))
    stop("These control parameters are not recognized: ",
         paste0(names(control)[is.na(control_match)], collapse = "\t"))

  control_default[control_match] <- control
  control <- control_default

  if(!control$criteria %in% c("delta_coef", "delta_likeli"))
    stop("Convergence criteria ", control$criteria, " not implemented")

  if(is.null(control$Q_0_term_for_fixed_E_step)){
    control$Q_0_term_for_fixed_E_step <- ifelse(
      control$method %in% c("UKF", "GMA"), 1, 1e5) # quite abritary values (especially the former - the latter is not too important)
  }

  if(!control$fixed_terms_method %in% c("M_step", "E_step"))
    stop("fixed_terms_method method '", control$fixed_terms_method,"' is not implemented")

  if(control$denom_term <= 0){
    stop("Method not implemented with penalty term (control$denom_term) equal to ", control$denom_term)
  } else if(control$denom_term < 1e-6)
    warning("Method is no tested with penalty term (control$denom_term) less than ", 1e-6)

  if(verbose)
    message("Finding Risk set")
  risk_set <-
    get_risk_obj(Y = X_Y$Y, by = by, max_T = ifelse(missing(max_T), max(X_Y$Y[X_Y$Y[, 3] == 1, 2]), max_T),
                 id = id, is_for_discrete_model = is_for_discrete_model)

  n_fixed <- ncol(X_Y$fixed_terms)
  est_fixed_in_E <- control$fixed_terms_method == "E_step" && n_fixed > 0

  if(missing(Q_0)){
    Q_0 = diag(10, n_params * order) # something large. Though depends on model, estimation method and data

    if(missing(Q))
      Q = diag(1, n_params) # (Very) arbitrary default
  }

  F_ <- get_F(order, n_params, n_fixed, est_fixed_in_E)

  # Check if there are any fixed coefficients. If not set the fixed
  # coefficients to an empty vector
  if(is.null(attr(X_Y$formula, "specials")[ddfixed_specials]))
    if(is.null(control$fixed_parems_start))
      control$fixed_parems_start <- vector("double")

  if(n_params == 0){
    # Model is fitted using ddhazard_fit_cpp for testing
    warning("The model can be estimated more effeciently by using get_survival_case_Weights_and_data and static_glm when there is no time varying parameters")
    a_0 = vector()

    if(is.null(control$fixed_parems_start))
      control$fixed_parems_start <- rep(0, ncol(X_Y$fixed_terms)) else
        control$fixed_parems_start <- control$fixed_parems_start

  } else if((missing_a_0 <- missing(a_0)) |
            (missing_fixed <- (is.null(control$fixed_parems_start)))){
    if(getOption("ddhazard_use_speedglm")){
      glm_func <- function(fam)
        suppressWarnings( # Get warning due to convergence failures when maxit = 1
          static_glm(formula = formula, data = data, max_T = max_T, risk_obj = risk_set,
                     maxit = 1, family = fam, speedglm = T))
    } else {
      glm_func <- function(fam)
        static_glm(formula = formula, data = data, max_T = max_T, risk_obj = risk_set,
                   control = stats::glm.control(epsilon = Inf), family = fam,
                   speedglm = F)
    }

    if(model == "logit"){
      tmp_mod = glm_func("binomial")

    } else if(model %in% exp_model_names){
      tmp_mod = glm_func("exponential")

    } else
      stop("Method not implemented to find initial values for '", model, "'. Please, provide intial values for a_0")

    is_fixed <-
      names(tmp_mod$coefficients) %in% colnames(X_Y$fixed_terms) |
      grepl("^ddFixed_intercept\\(", names(tmp_mod$coefficients), perl = TRUE)

    if(missing_a_0){
      message("a_0 not supplied. One iteration IWLS of static glm model is used")
      a_0 = rep(tmp_mod$coefficients[!is_fixed], order)
    }

    if(missing_fixed){
      control$fixed_parems_start <- tmp_mod$coefficients[is_fixed]
    }

    rm(tmp_mod)
  }

  if(is.vector(Q) && length(Q) == 1)
    Q <- matrix(Q)

  if(is.vector(Q_0) && length(Q_0) == 1)
    Q_0 <- matrix(Q_0)

  if(ncol(F_) != n_params * order + n_fixed * est_fixed_in_E)
    stop("F_ does not have the correct dimension. Its dimension should be ", n_params * order + n_fixed * est_fixed_in_E,
         " but it has ", ncol(F_), " columns")

  if(ncol(Q) != n_params)
    stop("Q does not have the correct dimension. Its dimension should be ", n_params,
         " but it has ", ncol(Q), " columns")

  if(ncol(Q_0) != n_params * order)
    stop("Q_0 does not have the correct dimension. Its dimension should be ", n_params * order,
         " but it has ", ncol(Q_0), " columns")

  if(length(a_0) != n_params * order)
    stop("a_0 does not have the correct length. Its length should be ", n_params * order,
         " but it has length ", length(a_0), " ")

  if(order > 1){
    tmp <- matrix(0., nrow = order * n_params, ncol = order * n_params)
    tmp[1:n_params, 1:n_params] <- Q
    Q <- tmp
  }

  if(est_fixed_in_E){
    # We need to add entries to the various matrices and vectors
    indicies_fix <- 1:n_fixed + n_params

    Q_new <- F_ # F_ already has the right dimensions
    Q_new[, ] <- 0
    Q_new[-indicies_fix, -indicies_fix] <- Q
    Q <- Q_new

    Q_0_new <- F_
    Q_0_new[, ] <- 0
    Q_0_new[-indicies_fix, -indicies_fix] <- Q_0
    if(length(indicies_fix) == 1)
      Q_0_new[indicies_fix, indicies_fix] <- control$Q_0_term_for_fixed_E_step else
        diag(Q_0_new[indicies_fix, indicies_fix]) <- control$Q_0_term_for_fixed_E_step
    Q_0 <- Q_0_new

    if(order == 1){
      a_0 <- c(a_0, control$fixed_parems_start)
    } else if(order == 2){
      first_half <- 1:(length(a_0)/2)
      a_0 <- c(a_0[first_half], control$fixed_parems_start, a_0[-first_half])
    } else
      stop("Order of '", order, "' not supported here")
    control$fixed_parems_start <- vector()
  }

  # Report pre-liminary stats
  if(verbose){
    tmp_tbl = matrix(NA_real_, nrow = risk_set$d, ncol = 2)
    colnames(tmp_tbl) = c("Risk size", "Num events")

    # Find the number of event in each bin
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

    # Insert in to table to print and print
    tmp_tbl[, "Num events"] <- n_events
    tmp_tbl[, "Risk size"] <- unlist(lapply(risk_set$risk_sets, length))

    message("Total number of included events are ", sum(tmp_tbl[, 2]), " of ",   sum(X_Y$Y[, 3]))
    message("Size of risk set and number of events in each risk set are ([Risk set size]:[# events]):")

    tmp_tbl[, 1] <- sprintf("%8s", tmp_tbl[, 1])
    tmp_tbl[, 2] <- sprintf("%-8s", tmp_tbl[, 2])
    tmp_message <- sprintf("%21s", apply(tmp_tbl, 1, paste, collapse = " : "))
    msg_final <- tmp_message[1]
    for(i in seq_along(tmp_message)[-1])
      msg_final <- paste0(msg_final, if((i - 1) %% 4 > 0) " " else "\n", tmp_message[i])
    message(msg_final)

    message("Running EM")
  }

  if(control$permu)
    eval(get_permu_data_exp(X_Y[1:3], risk_set, weights))

  if(!est_fixed_in_E){
    X_Y$X <- t(X_Y$X) # we transpose for performance due to the column-major
    # ordering. The logic is that we primary look up indviduals
    # (i.e. columns in the tranpose)
    X_Y$fixed_terms <- t(X_Y$fixed_terms) # same
  } else {
    # We add the firm covariates to the design matrix. Later, we recover the fixed covariates
    X_Y$X <- t(cbind(X_Y$X, X_Y$fixed_terms))
    X_Y$fixed_terms <- matrix(nrow = 0, ncol = ncol(X_Y$X))
  }

  result <- NA
  k_vals <- seq_len(control$LR_max_try) - 1
  for(k in k_vals){
    tryCatch({
      result <- ddhazard_no_validation(a_0 = a_0, Q_0 = Q_0, F_ = F_, verbose = verbose, Q = Q,
                                       risk_set= risk_set, X_Y = X_Y, order = order, model = model,
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

  if(k != 0)
    message("Did not succed to fit the model wit a learning rate of ", control$LR,
            ". The learning rate was decrease by a factor ", control$LR_decrease_fac^k, " to yield a fit")

  if(model != "logit" || control$method != "EKF")
    result$lag_one_cov = NULL # only checked the value for the logit model with EKF

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

    result$lag_one_cov <- result$lag_one_cov[-indicies_fix, -indicies_fix, , drop = F]

    result$Q  <- result$Q[-indicies_fix, -indicies_fix, drop = F]

    F_ <- F_[-indicies_fix, -indicies_fix, drop = F]
  }

  if(control$permu)
    eval(get_permu_data_rev_exp(list(), risk_set, weights))

  # Set names
  tmp_names = rep(rownames(X_Y$X), order)
  colnames(result$a_t_d_s) = tmp_names
  dimnames(result$V_t_d_s) = list(tmp_names, tmp_names, NULL)
  dimnames(result$Q) = list(tmp_names, tmp_names)

  result$fixed_effects <- c(result$fixed_effects)
  names(result$fixed_effects) <- fixed_names <- gsub(
    "(^ddFixed\\()(.+)(\\)$)","\\2", rownames(X_Y$fixed_terms))

  if(est_fixed_in_E){
    dimnames(result$Q_0) = list(c(tmp_names, fixed_names), c(tmp_names, fixed_names))
  } else
    dimnames(result$Q_0) = dimnames(result$Q)

  if(model == "logit") {
    res <- list(
      hazard_func =  function(eta, ...){
        exp_ = exp(eta)
        exp_/(1 + exp_)
      },
      hazard_first_deriv = function(beta, x_, ...){
        exp_ = exp(beta %*% x_)
        x_ * exp_ / (exp_ + 1)^2
      },
      var_func = function(eta, ...){
        exp_ = exp(eta)
        exp_ / (1 + exp_)^2
      })

  }else if(model %in% exp_model_names){
    res <- list(
    hazard_func =  function(eta, tstart, tstop, ...){
      1 - exp( - exp(eta) * (tstop - tstart))
    },
    hazard_first_deriv = function(beta, x_, tstart, tstop, ...){
      eta <- beta %*% x_
      x_ * (tstop - tstart) * exp(eta - exp(eta) * (tstop - tstart))
    },
    var_func = NULL)
  }

  structure(c(
    res, list(
    formula = formula,
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
    id = if(control$save_data) id else NULL,
    order = order, F_ = F_,
    method = control$method,
    model = model,
    est_Q_0 = control$est_Q_0,
    control = control),
    LR = LR),
    "class" = "fahrmeier_94")
}

ddhazard_no_validation <- function(a_0, Q_0, F_, verbose, Q,
                                   risk_set, X_Y, order, model, LR,
                                   n_fixed_terms_in_state_vec,
                                   weights = weights,
                                   control){
  ddhazard_fit_cpp(a_0 = a_0, Q_0 = Q_0, F_ = F_, verbose = verbose,
                   Q = Q, n_max = control$n_max,
                   risk_obj = risk_set, eps = control$eps,
                   X = X_Y$X, fixed_terms = X_Y$fixed_terms,
                   fixed_parems_start = control$fixed_parems_start,
                   tstart = X_Y$Y[, 1], tstop = X_Y$Y[, 2],
                   order_ = order,
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
                   EKF_batch_size = control$EKF_batch_size)
}


get_F <- function(order, n_params, n_fixed, est_fixed_in_E){
  if(order == 1){
    return(diag(rep(1, n_params + n_fixed * est_fixed_in_E)))
  }
  else if(order == 2){
    indicies_cur <- 1:n_params
    indicies_lag <- n_params + 1:n_params + ifelse(est_fixed_in_E, n_fixed, 0)
    indicies_fix <- if(est_fixed_in_E) 1:n_fixed + n_params else vector()

    F_ = matrix(NA_real_,
                nrow = 2 * n_params + n_fixed * est_fixed_in_E,
                ncol = 2 * n_params + n_fixed * est_fixed_in_E)
    F_[indicies_cur, indicies_cur] = diag(2, n_params)
    F_[indicies_lag, indicies_cur] = diag(1, n_params)
    F_[indicies_cur, indicies_lag] = diag(-1, n_params)
    F_[indicies_lag, indicies_lag] = 0

    if(length(indicies_fix) > 0){
      F_[indicies_fix, ] <- 0
      F_[, indicies_fix] <- 0
      if(length(indicies_fix) > 1)
        diag(F_[indicies_fix, indicies_fix]) <- 1 else
          F_[indicies_fix, indicies_fix] <- 1
    }

    return(F_)
  } else stop("Method not implemented for order ", order)
}

exp_model_names <- c("exp_bin",
                     "exp_clip_time", "exp_clip_time_w_jump")
