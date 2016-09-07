#' Residuals function for result of ddhazard fit
#' @export
residuals.fahrmeier_94 = function(object, type = c("stdSpaceErrors", "Pearson"), data_){

  if(!object$model %in% c("logit"))
    stop("Functions for model '",  object$model, "' is not implemented")

  type = type[1]

  if(type == "stdSpaceErrors"){
    warning("Std state space error is not tested")

    if(!object$method %in% c("EKF"))
      stop("Functions for stdSpaceErrors with method '", object$method, "' is not implemented")

    if(object$order != 1)
      stop("Method for stdSpaceErrors is not implemented for order ", object$order)

    res_ = diff(object$a_t_d_s)
    vars_ = object$V_t_d_s
    G = array(NA_real_, dim = dim(object$lag_one_cor))
    for(i in seq_len(dim(vars_)[3] - 1))
      G[, , i] = chol(
        vars_[, , i + 1] + vars_[, , i] - object$lag_one_cor[, , i] - t(object$lag_one_cor[, , i]))

    res_std = matrix(NA_real_, nrow = dim(res_)[1], ncol = dim(res_)[2])
    for(i in seq_len(dim(vars_)[3] - 1))
      res_std[i, ] = res_[i, ] %*% solve(G[, , i])

    return(structure(res_std,
                     "class" = "fahrmeier_94_stdSpaceErrors"))
  }
  if(type == "Pearson" || type == "Raw"){
    if(missing(data_) || is.na(object$risk_set))
      stop("Missing risk set or data")

    # Wee need these to check if there is an event in the bin
    new_stop = object$risk_set$stop_new
    new_event = object$risk_set$new_events_flags
    res = list()

    for(i in seq_along(object$risk_set$risk_sets)){
      start_ = object$times[i]
      stop_ = object$times[i + 1]
      r_set = object$risk_set[[1]][[i]]

      tmp_dat = data_[r_set, ]
      tmp_dat$tstart = rep(start_, length(r_set))
      tmp_dat$tstop = rep(stop_, length(r_set))
      p_est = predict(object, tmp_dat, tstart = "tstart", tstop = "tstop")$fits

      Y = (new_stop[r_set] == stop_) * new_event[r_set]

      if(type == "Raw"){
        res[[i]] = cbind(Raw_res = Y - p_est,
                         p_est = p_est,
                         Y = Y,
                         row_num = r_set)
        next
      }

      res[[i]] = cbind(Pearson_res = (Y - p_est)/sqrt(p_est * (1 - p_est)),
                       p_est = p_est,
                       Y = Y,
                       row_num = r_set)
    }

    if(type == "Raw")
      return(structure(res,
                       "class" = "fahrmeier_94_Raw_res"))

    return(structure(res,
                     "class" = "fahrmeier_94_Pearson"))
  }
}
