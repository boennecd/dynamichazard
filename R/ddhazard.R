#' Function to fit dynamic discrete hazard model
#' @param formula \code{coxph} like formula with \code{survival::Surv(tstart, tstop, event)} on the left hand site of \code{~}
#' @parem data Data frame or enviroment containing the outcome and co-variates
#' @parem model "logit" or "exponential" for the discrete time function using the logistic function or for the continous time model with exponential arrival times
#' @param by Interval length of the bins in which parameters are fixed
#' @parem max_T End of the last interval. The last stop time with an event is selected if the parameter is omitted
#' @param id Vector of ids for each row of the in the design matrix
#' @parem a_0 Vector with \eqn{\vec{a}_0} for the first iteration (optional)
#' @param Q_0 Covariance matrix for the prior \eqn{\mathbf{Q}_0}
#' @parem Q Initial covariance matrix \eqn{\mathbf{Q}} for the state equation
#' @parem order_ Order of the random walk
#' @parem control list of control variables
#' @export
ddhazard = function(formula, data,
                    model = "logit",
                    by, max_T, id,
                    a_0, Q_0, Q = Q_0,
                    order_ = 1, control = list(),
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
    Q_0 = diag(10, n_parems * order_) # something large. TODO: what is large?

    if(missing(Q))
      Q = diag(c(rep(1, n_parems), rep(0, n_parems * (order_ - 1)))) # TODO: What to set here?
  }

  if(order_ == 1){
    F_ = diag(1, n_parems)
  }
  else if(order_ == 2){
    F_ = matrix(NA_real_, nrow = 2 * n_parems, ncol = 2 * n_parems)
    F_[1:n_parems, 1:n_parems] = diag(2, n_parems)
    F_[n_parems + 1:n_parems, 1:n_parems] = diag(1, n_parems)
    F_[1:n_parems, n_parems + 1:n_parems] = diag(-1, n_parems)
    F_[n_parems + 1:n_parems, n_parems + 1:n_parems] = 0
  } else stop("Method not implemented for order ", order_)

  if(model == "logit"){
    is_for_discrete_model <- TRUE
  } else if (model == "exponential"){
    is_for_discrete_model <- FALSE
  } else
    stop("Model '", model, "' is not implemented")

  control_default <- list(kappa = NULL, alpha = NULL, beta = NULL,
                          NR_eps = NULL, LR = NULL, n_max = 10^2, eps = 10^-3,
                          est_Q_0 = F, method = "EKF", save_risk_set = T)
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

  if(missing(a_0) && model == "logit"){
    # Assume that logit models is used
    message("a_0 not supplied. One iteration IWLS of static logit model is used")
    tmp_mod = static_glm(form = formula, data = data, risk_obj = risk_set,
                         control = glm.control(epsilon = Inf), family = "binomial")
    a_0 = rep(tmp_mod$coefficients, order_)
    rm(tmp_mod)

  } else if (missing(a_0) && model == "exponential"){
    message("a_0 not supplied. One iteration IWLS of static glm model is used")
    tmp_mod = static_glm(form = formula, data = data, max_T = max_T,
                         control = glm.control(epsilon = Inf), family = "exponential")
    a_0 = rep(tmp_mod$coefficients, order_)
    rm(tmp_mod)
  }

  if(ncol(F_) != n_parems * order_ ||
     ncol(Q) != n_parems * order_ ||
     ncol(Q_0) != n_parems * order_ ||
     length(a_0) != n_parems * order_)
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

  result = ddhazard_fit_cpp_prelim(a_0 = a_0, Q_0 = Q_0, F_ = F_, verbose = verbose,
                                   Q = Q, n_max = control$n_max,
                                   risk_obj = risk_set, eps = control$eps, X = X_Y$X,
                                   tstart = X_Y$Y[, 1], tstop = X_Y$Y[, 2],
                                   order_ = order_,
                                   est_Q_0 = control$est_Q_0, method = control$method,
                                   model = model,
                                   kappa = control$kappa, alpha = control$alpha, beta = control$beta,
                                   NR_eps = control$NR_eps,
                                   LR = control$LR)

  # Set names
  tmp_names = rep(colnames(X_Y$X), order_)
  colnames(result$a_t_d_s) = tmp_names
  dimnames(result$V_t_d_s) = list(tmp_names, tmp_names, NULL)
  dimnames(result$Q) = dimnames(result$Q_0) = list(tmp_names, tmp_names)

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
    hazard_first_deriv = function(beta, x_, ...){
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
    n_iter = result$n_iter,
    Q = result$Q,
    Q_0 = result$Q_0,
    n_risk = unlist(lapply(risk_set$risk_sets, length)),
    times = c(min(X_Y$Y[, 1]), risk_set$event_times),
    risk_set = if(control$save_risk_set) risk_set else NA,
    order = order_, F_ = F_,
    method = control$method,
    model = model,
    est_Q_0 = control$est_Q_0)),
    "class" = "fahrmeier_94")
}
