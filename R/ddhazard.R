#' Function to fit dynamic discrete hazard model
#' @param formula \code{\link[survival]{coxph}} like formula with \code{\link[survival]{Surv}(tstart, tstop, event)} on the left hand site of \code{~}
#' @param data Data frame or enviroment containing the outcome and co-variates
#' @param model "logit" or "exponential" for the discrete time function using the logistic function or for the continous time model with exponential arrival times
#' @param by Interval length of the bins in which parameters are fixed
#' @param max_T End of the last interval. The last stop time with an event is selected if the parameter is omitted
#' @param id Vector of ids for each row of the in the design matrix
#' @param a_0 Vector with \eqn{\vec{a}_0} for the first iteration (optional). Default is estimates from static model (see \code{\link{static_glm}})
#' @param Q_0 Covariance matrix for the prior \eqn{\mathbf{Q}_0}
#' @param Q Initial covariance matrix \eqn{\mathbf{Q}} for the state equation
#' @param order Order of the random walk
#' @param control List of control variables (see details below)
#' @param verbose \code{TRUE} if you want status messages during execution
#'
#' @details
#' Function to estimate a binary regression function of the form where the regression parameters follows a given order random walk. The order is specified by the \code{order} argument. 1. and 2. order random walks is implemented. The regression parameters are updated at time \code{by}, 2\code{by}, ..., \code{max_T}. See the vignette 'ddhazard' for more details
#'
#' The Extended Kalman filter or Uncented Kalman filter needs an initial co-variance matrix \code{Q_0} and state vector \code{a_0}. An estimate from a time-invariant model is provided for \code{a_0} if it is not supplied (the same model you would get from \code{\link{static_glm}} function). A diagonal matrix with large entries is recommended for \code{Q_0}. What is large dependents on the data set and \code{model}. Further, a variance matrix for the first iteration \code{Q} is needed. It is recommended to select diagonal matrix with low values for the latter. The \code{Q}, \code{a_0} and optionally \code{Q_0} is estimated with an EM-algorithm
#'
#' The model is specified through the \code{model} argument. Currently, \code{'logit'} and \code{'exponential'} is provided. The former uses an logistic model where outcomes are binned into the intervals. Be aware that there can be loss of information due to the binning. It is key for the logit model that the \code{id} argument is provided if individuals in the data set have time varying co-variates. The latter model use an exponential model for the arrival times where there is no loose information due to binning
#'
#' @section Control:
#' The \code{control} argument allows you to pass a \code{list} to select additional parameters. See the vignette 'ddhazard' for more information on hyper parameters, \code{LR} and \code{NR_eps}
#' \tabular{ll}{
#' \code{method}\verb{  }\tab Set to the method to use in the E-step. Either \code{"EKF"} for the Extended Kalman Filter or \code{"UKF"}for the Unscented Kalman Filter. \code{"EKF"} is the default \cr
#' \code{LR}\verb{  }\tab Learning rate for the Extended Kalman filter. Default is 1 \cr
#' \code{NR_eps}\verb{  }\tab Tolerance for the Extended Kalman filter. Default is \code{NULL} which means that no extra iteration is made in the correction step \cr
#' \code{alpha}\verb{  }\tab Hyper parameter \eqn{\alpha} in the Unscented Kalman Filter. Default is 1 \cr
#' \code{beta}\verb{  }\tab Hyper parameter \eqn{\beta} in the Unscented Kalman Filter. Default is 2 \cr
#' \code{kappa}\verb{  }\tab Hyper parameter \eqn{\kappa} in the Unscented Kalman Filter. Default is 0 \cr
#' \code{n_max}\verb{  }\tab Maximum number of iteration in the EM-algorithm \cr
#' \code{eps}\verb{  }\tab Tolerance parameter for the EM-algorithm\cr
#' \code{est_Q_0}\verb{  }\tab\code{TRUE} if you want the EM-algorithm to estimate \code{Q_0}. Default is \code{FALSE}\cr
#' \code{save_risk_set}\verb{  }\tab \code{TRUE} if you want to save the list from \code{\link{get_risk_obj}} used to estimate the model. It may be needed for later call to \code{residuals}, \code{plot} and \code{logLike}. Can be set to \code{FALSE} to save memory\cr
#' \code{save_data}\verb{  }\tab \code{TRUE} if you want to save the list \code{data} argument. It may be needed for later call to \code{residuals}, \code{plot} and \code{logLike}. Can be set to \code{FALSE} to save memory
#'}
#'
#' @return
#' A list with class \code{fahrmeier_94}. The list contains:
#' \tabular{ll}{
#' \code{formula}\verb{ }\tab The passed formula \cr
#' \code{state_vecs}\verb{ }\tab 2D matrix with the estimated state vectors (regression parameters) in each bin \cr
#' \code{state_vars}\verb{ }\tab 3D matrix with smoothed variance estimates for each state vector \cr
#' \code{lag_one_cor}\verb{ }\tab 3D Matrix with lagged correlation matrix for each for each change in the state vector \cr
#' \code{n_risk}\verb{ }\tab The number of observations in each interval \cr
#' \code{times}\verb{ }\tab The interval borders \cr
#' \code{risk_set}\verb{ }\tab The object from \code{\link{get_risk_obj}} if saved \cr
#' \code{data}\verb{ }\tab The \code{data} argument if saved \cr
#' \code{order}\verb{ }\tab Order of the random walk \cr
#' \code{F_}\verb{ }\tab Matrix with that map transition from one state vector to the next \cr
#' \code{method}\verb{ }\tab Method used in the E-step \cr
#' \code{est_Q_0}\verb{ }\tab\code{TRUE} if \code{Q_0} was estimated in the EM-algorithm \cr
#' \code{hazard_func}\verb{ }\tab Hazard function \cr
#' \code{hazard_first_deriv}\verb{ }\tab First derivative of the hazard function with respect to the linear predictor \cr
#'}
#'
#' @seealso
#' \code{\link[=plot.fahrmeier_94]{plot}}, \code{\link[=residuals.fahrmeier_94]{residuals}}, \code{\link[=predict.fahrmeier_94]{predict}}, \code{\link{static_glm}}
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
                    order = 1, control = list(),
                    verbose = F){
  X_Y = get_design_matrix(formula, data)
  n_parems = ncol(X_Y$X)
  tmp_n_failures = sum(X_Y$Y[, 3])

  if(missing(id)){
    if(verbose)
      warning("You did not parse and ID argument. I do not hink this is what you want ...")
    id = 1:nrow(data)
  }

  if(missing(Q_0)){
    Q_0 = diag(10, n_parems * order) # something large. Though depends on model, estimation method and data

    if(missing(Q))
      Q = diag(c(rep(1, n_parems), rep(0, n_parems * (order - 1)))) # (Very) arbitrary default
  }

  if(order == 1){
    F_ = diag(1, n_parems)
  }
  else if(order == 2){
    F_ = matrix(NA_real_, nrow = 2 * n_parems, ncol = 2 * n_parems)
    F_[1:n_parems, 1:n_parems] = diag(2, n_parems)
    F_[n_parems + 1:n_parems, 1:n_parems] = diag(1, n_parems)
    F_[1:n_parems, n_parems + 1:n_parems] = diag(-1, n_parems)
    F_[n_parems + 1:n_parems, n_parems + 1:n_parems] = 0

  } else stop("Method not implemented for order ", order)

  if(model == "logit"){
    is_for_discrete_model <- TRUE

  } else if (model == "exponential"){
    is_for_discrete_model <- FALSE

  } else
    stop("Model '", model, "' is not implemented")

  control_default <- list(kappa = NULL, alpha = NULL, beta = NULL,
                          NR_eps = NULL, LR = 1, n_max = 10^2, eps = 10^-3,
                          est_Q_0 = F, method = "EKF", save_risk_set = T,
                          save_data = T, eps_fixed_parems = 1e-3,
                          max_it_fixed_parems = 10, fixed_effect_chunk_size = 1e4,
                          debug = F, fixed_parems_start = NULL, LR_max_try = 10,
                          LR_decrease_fac = 1.5, EKF_inv_Cov_method = "org")

  if(any(is.na(control_match <- match(names(control), names(control_default)))))
    stop("These control parameters are not recognized: ",
         paste0(names(control)[is.na(control_match)], collapse = "\t"))

  control_default[control_match] <- control
  control <- control_default

  if(verbose)
    message("Finding Risk set")
  risk_set <-
    get_risk_obj(Y = X_Y$Y, by = by, max_T = ifelse(missing(max_T), max(X_Y$Y[X_Y$Y[, 3] == 1, 2]), max_T),
                 id = id, is_for_discrete_model = is_for_discrete_model)

  if(n_parems == 0){
    # Model is fitted using ddhazard_fit_cpp for testing
    warning("The model can be estimated more effeciently by using get_survival_case_Weigths_and_data and static_glm when there is no time varying parameters")
    a_0 = vector()

    if(is.null(control$fixed_parems_start))
      control$fixed_parems_start <- rep(0, ncol(X_Y$fixed_terms)) else
        control$fixed_parems_start <- control$fixed_parems_start

  } else if((missing_a_0 <- missing(a_0)) |
            (missing_fixed <- is.null(control$fixed_parems_start))){
    if(model == "logit"){
      tmp_mod = static_glm(formula = formula, data = data, risk_obj = risk_set,
                           control = stats::glm.control(epsilon = Inf), family = "binomial")

    } else if(model == "exponential"){
      tmp_mod = static_glm(formula = formula, data = data, max_T = max_T,
                           control = stats::glm.control(epsilon = Inf), family = "exponential")

    } else
      stop("Method not implemented to find initial values for '", model, "'")

    is_fixed <- seq_along(tmp_mod$coefficients) %in% (attr(X_Y$formula, "specials")$ddFixed - 1)

    if(missing_a_0){
      message("a_0 not supplied. One iteration IWLS of static glm model is used")
      a_0 = rep(tmp_mod$coefficients[!is_fixed], order)

    }
    if(missing_fixed){
      control$fixed_parems_start <- tmp_mod$coefficients[is_fixed]

    }

    rm(tmp_mod)
  }

  if(ncol(F_) != n_parems * order ||
     ncol(Q) != n_parems * order ||
     ncol(Q_0) != n_parems * order ||
     length(a_0) != n_parems * order)
    stop("One of the input vector or matrices do not match with the order and number of parameters")

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

    message("Total number of included events are ", sum(tmp_tbl[, 2]), " of ", tmp_n_failures)
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

  X_Y$X <- t(X_Y$X) # we transpose for performance due to the column-major
                    # ordering. The logic is that we primary look up indviduals
                    # (i.e. columns in the tranpose)
  X_Y$fixed_terms <- t(X_Y$fixed_terms) # same

  result <- NA
  k_vals <- if(control$method == "EKF") seq_len(control$LR_max_try) - 1 else 0
  for(k in k_vals){
    tryCatch({
      result <- ddhazard_fit_cpp(a_0 = a_0, Q_0 = Q_0, F_ = F_, verbose = verbose,
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
                                 LR = control$LR / control$LR_decrease_fac^(k),
                                 eps_fixed_parems = control$eps_fixed_parems,
                                 max_it_fixed_parems = control$max_it_fixed_parems,
                                 fixed_effect_chunk_size = control$fixed_effect_chunk_size,
                                 debug = control$debug,
                                 EKF_inv_Cov_method = control$EKF_inv_Cov_method)
    }, error = function(e)
      if(!grepl("^ddhazard_fit_cpp estimation error:", e$message))
        stop(e))

    LR <- control$LR / control$LR_decrease_fac^(k)
    if(did_fit <- !(length(result) == 1 && is.na(result)))
      break
  }

  # Check if fit succeeded
  if(!did_fit)
    stop("Failed to estimate model. The following learning rates was tried: ",
         paste0(control$LR / control$LR_decrease_fac^k_vals, collapse = ", "),
         ". Try decreasing the learning rate")

  if(k != 0)
    message("Did not succed to fit the model wit a learning rate of ", control$LR,
            ". The learning rate was decrease by a factor ", control$LR_decrease_fac^k, " to yield a fit")

  # Set names
  tmp_names = rep(rownames(X_Y$X), order)
  colnames(result$a_t_d_s) = tmp_names
  dimnames(result$V_t_d_s) = list(tmp_names, tmp_names, NULL)
  dimnames(result$Q) = dimnames(result$Q_0) = list(tmp_names, tmp_names)

  result$fixed_effects <- c(result$fixed_effects)
  names(result$fixed_effects) <-gsub("(^ddFixed\\()(.+)(\\)$)","\\2", rownames(X_Y$fixed_terms))

  if(model == "logit") {
    res <- list(
    hazard_func =  function(eta, ...){
      exp_ = exp(eta)
      exp_/(1 + exp_)
      },
      hazard_first_deriv = function(beta, x_, ...){
        exp_ = exp(beta %*% x_)
        x_ * exp_ / (exp_ + 1)^2
      })

  }else if(model == "exponential"){
    res <- list(
    hazard_func =  function(eta, tstart, tstop, ...){
      1 - exp( - exp(eta) * (tstop - tstart))
    },
    hazard_first_deriv = function(beta, x_, tstart, tstop, ...){
      eta <- beta %*% x_
      x_ * (tstop - tstart) * exp(eta - exp(eta) * (tstop - tstart))
    })
  }

  structure(c(
    res, list(
    formula = X_Y$formula,
    state_vecs = result$a_t_d_s,
    state_vars = result$V_t_d_s,
    lag_one_cor = result$lag_one_cor,
    fixed_effects = result$fixed_effects,
    n_iter = result$n_iter,
    Q = result$Q,
    Q_0 = result$Q_0,
    n_risk = unlist(lapply(risk_set$risk_sets, length)),
    times = c(min(X_Y$Y[, 1]), risk_set$event_times),
    risk_set = if(control$save_risk_set) risk_set else NULL,
    data = if(control$save_data) data else NULL,
    order = order, F_ = F_,
    method = control$method,
    model = model,
    est_Q_0 = control$est_Q_0),
    LR = LR),
    "class" = "fahrmeier_94")
}
