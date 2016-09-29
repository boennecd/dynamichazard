#' Residuals function for result of ddhazard fit
#' @export
residuals.fahrmeier_94 = function(object, type = c("std_space_error", "space_error", "pearson", "raw"), data_){
  if(!object$model %in% c("logit"))
    stop("Functions for model '",  object$model, "' is not implemented")

  type = type[1]

  if(type %in% c("std_space_error", "space_error")){
    return(space_errors(object, data_, type == "std_space_error"))
  }

  if(type == "pearson" || type == "raw"){
    return(obs_res())
  }

  stop("Method '", type, "' not implemented for residuals method")
}

space_errors <- function(object, data_, standardize){
  warning("Std state space error is not tested")

  if(!object$method %in% c("EKF"))
    stop("Functions for with method '", object$method, "' is not implemented")

  if(object$order != 1)
    stop("Method for is not implemented for order ", object$order)

  res = diff(object$a_t_d_s)

  covs <- vars <- object$V_t_d_s
  for(i in seq_len(dim(vars)[3] - 1))
    covs[, , i] = vars[, , i + 1] + vars[, , i] - object$lag_one_cor[, , i] - t(object$lag_one_cor[, , i])

  if(standardize){
    for(i in seq_len(dim(vars)[3] - 1))
      res[i, ] = res[i, ] %*% solve(chol(covs[, , i]))
  }

  return(structure(list(residuals = res,
                        standardize = standardize,
                        Covariances = covs),
                   "class" = "fahrmeier_94_SpaceErrors"))
}

obs_res <- function(...){
  # if(missing(data_) || is.na(object$risk_set))
  #   stop("Missing risk set or data")
  #
  # # Wee need these to check if there is an event in the bin
  # new_stop = object$risk_set$stop_new
  # new_event = object$risk_set$new_events_flags
  # res = list()
  #
  # for(i in seq_along(object$risk_set$risk_sets)){
  #   start_ = object$times[i]
  #   stop_ = object$times[i + 1]
  #   r_set = object$risk_set[[1]][[i]]
  #
  #   tmp_dat = data_[r_set, ]
  #   tmp_dat$tstart = rep(start_, length(r_set))
  #   tmp_dat$tstop = rep(stop_, length(r_set))
  #   p_est = predict(object, tmp_dat, tstart = "tstart", tstop = "tstop")$fits
  #
  #   Y = (new_stop[r_set] == stop_) * new_event[r_set]
  #
  #   if(type == "Raw"){
  #     res[[i]] = cbind(Raw_res = Y - p_est,
  #                      p_est = p_est,
  #                      Y = Y,
  #                      row_num = r_set)
  #     next
  #   }
  #
  #   res[[i]] = cbind(Pearson_res = (Y - p_est)/sqrt(p_est * (1 - p_est)),
  #                    p_est = p_est,
  #                    Y = Y,
  #                    row_num = r_set)
  # }
  #
  # if(type == "Raw")
  #   return(structure(res,
  #                    "class" = "fahrmeier_94_Raw_res"))
  #
  # return(structure(res,
  #                  "class" = "fahrmeier_94_Pearson"))
}
