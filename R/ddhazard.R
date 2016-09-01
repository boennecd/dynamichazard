#' Function to fit dynamic discrete hazard model
#' @export
ddhazard = function(formula, data, by,
                    a_0,
                    Q_0,
                    id,
                    Q = Q_0,
                    n_max = 10^2, eps = 10^-3,
                    save_all_output = F,
                    max_T, save_risk_set = F,
                    order_ = 1,
                    est_Q_0 = T,
                    verbose = F){
  X_Y = benssurvutils::get_design_matrix(formula, data)
  n_parems = ncol(X_Y$X)
  tmp_n_failures = sum(X_Y$Y[, 3])

  if(missing(a_0)){
    # Assume that logit models is used
    message("a_0 not supplied. One iteration IWLS of static logit model is used")
    tmp_mod = glm(X_Y$Y[, 3] ~ X_Y$X - 1, # design mat already have intercept
                  family = binomial, control = glm.control(epsilon = Inf))
    a_0 = rep(tmp_mod$coefficients, order_)
    rm(tmp_mod)
  }

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
    # TODO: Test
    F_ = matrix(NA_real_, nrow = 2 * n_parems, ncol = 2 * n_parems)
    F_[1:n_parems, 1:n_parems] = diag(2, n_parems)
    F_[n_parems + 1:n_parems, 1:n_parems] = diag(1, n_parems)
    F_[1:n_parems, n_parems + 1:n_parems] = diag(-1, n_parems)
    F_[n_parems + 1:n_parems, n_parems + 1:n_parems] = 0
  } else stop("Method not implemented for order ", order_)

  if(ncol(F_) != n_parems * order_ ||
     ncol(Q) != n_parems * order_ ||
     ncol(Q_0) != n_parems * order_ ||
     length(a_0) != n_parems * order_)
    stop("One of the input vector or matrices do not match with the order and number of parameters")

  if(verbose)
    message("Finding Risk set")
  risk_set <-
    benssurvutils::get_risk_sets(Y = X_Y$Y, by = by, max_T = ifelse(missing(max_T), max(X_Y$Y[X_Y$Y[, 3] == 1, 2]), max_T), id = id)

  X_Y$Y[, 2] = risk_set$stop_new # Update to new rounded risk stop times
  X_Y$Y[, 3] = risk_set$new_events_flags # update event flags

  # Report pre-liminary stats
  if(verbose){
    tmp_tbl = matrix(NA_real_, nrow = risk_set$d, ncol = 2)
    colnames(tmp_tbl) = c("Risk size", "Num events")
    tmp_stop = min(X_Y$Y[, 1])
    for(i in seq_len(risk_set$d)){
      tmp_stop = tmp_stop + risk_set$I_len[i]
      r_set = risk_set$risk_sets[[i]]
      tmp_tbl[i, ] = c(length(r_set), sum(X_Y$Y[r_set, 3] * (X_Y$Y[r_set, 2] == tmp_stop)))
    }
    message("Size of risk set and number of events in each risk set are:")
    message(paste(apply(tmp_tbl, 1, paste, collapse = " : "), collapse = ",\t"))
    message("Total number of included events are ", sum(tmp_tbl[, 2]), " of ", tmp_n_failures)
    rm(tmp_tbl, tmp_stop, tmp_n_failures)

    message("Running EM")
  }

  result = ddhazard_fit_cpp_prelim(a_0 = a_0, Q_0 = Q_0, F_ = F_, verbose = verbose, save_all_output = save_all_output,
                                   Q = Q, n_max = n_max,
                                   risk_obj = risk_set, eps = eps, X = X_Y$X,
                                   tstart = X_Y$Y[, 1], tstop = X_Y$Y[, 2], events = X_Y$Y[, 3],
                                   order_ = order_,
                                   est_Q_0 = est_Q_0)

  # Set names
  tmp_names = rep(colnames(X_Y$X), order_)
  colnames(result$a_t_d_s) = tmp_names
  dimnames(result$V_t_d_s) = list(tmp_names, tmp_names, NULL)
  dimnames(result$Q) = dimnames(result$Q_0) = list(tmp_names, tmp_names)

  structure(list(
    formula = X_Y$formula,
    a_t_d_s = result$a_t_d_s,
    V_t_d_s = result$V_t_d_s,
    lag_one_cor = result$lag_one_cor,
    n_iter = result$n_iter,
    Q = result$Q,
    Q_0 = result$Q_0,
    n_risk = unlist(lapply(risk_set$risk_sets, length)),
    times = c(min(X_Y$Y[, 1]), risk_set$event_times),
    hazard_func =  function(eta){
      exp_ = exp(eta)
      exp_/(1 + exp_)
    },
    hazard_first_deriv = function(beta, x_){
      exp_ = exp(beta %*% x_)
      x_ * exp_ / (exp_ + 1)^2
    },
    risk_set = if(save_risk_set) risk_set else NA,
    order = order_, F_ = F_),
    "class" = "fahrmeier_94")
}
