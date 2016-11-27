#' Residuals function for the result of \code{\link{ddhazard}} fit
#'
#' @param object Result of \code{\link{ddhazard}} call
#' @param type Type of residuals. Four possible values: \code{"std_space_error"}, \code{"space_error"}, \code{"pearson"} and \code{"raw"}. See below for details
#' @param data Data frame with data for Pearson or raw residuals
#' @param ... Not used
#'
#' @section Pearson and raw residuals:
#' Is the result of a call with a \code{type} argument of either \code{"pearson"} or \code{"raw"} for Pearson residuals or raw residuals. Returns a list with class \code{"fahrmeier_94_res"} with the following elements
#' \tabular{ll}{
#' \code{residuals}\verb{ }\tab List of residuals for each bin. Each element of the list contains a 2D array where the rows corresponds to the passed \code{data} and columns are the residuals (\code{residuals}), estimated probability of death (\code{p_est}), outcome (\code{Y}) and row number in the initial dataset (\code{row_num}). The \code{data} rows will only have a residuals in a given risk list if they are at risk in that risk set\cr
#' \code{type}\verb{ }\tab The type of residual\cr
#'}
#'
#' @section State space errors:
#' Is the result of a call with a \code{type} argument of either \code{"std_space_error"} or \code{"space_error"}. Returns a list with class \code{"fahrmeier_94_SpaceErrors"} with the following elements
#' \tabular{ll}{
#' \code{residuals}\verb{ }\tab 2D array with either standardised or unstandardised state space errors. The row are bins and the columns are the parameters in the regression \cr
#' \code{standardize}\verb{ }\tab \code{TRUE} if standardised state space errors \cr
#' \code{Covariances}\verb{ }\tab 3D array with the smoothed co-variance matrix for each set of the state space errors \cr
#'}
#'
#' @export
residuals.fahrmeier_94 = function(object, type = c("std_space_error", "space_error", "pearson", "raw"), data = NULL, ...){
  type = type[1]

  if(type %in% c("std_space_error", "space_error")){
    return(space_errors(object, data, type == "std_space_error"))
  }

  if(!object$model %in% c("logit", "exponential"))
    stop("Functions for model '",  object$model, "' is not implemented")

  if(type == "pearson" || type == "raw"){
    return(obs_res(object, if(is.null(object$data)) data else object$data, type))
  }

  stop("Method '", type, "' not implemented for residuals method")
}

space_errors <- function(object, data, standardize){
  if(!object$method %in% c("EKF"))
    stop("Functions for with method '", object$method, "' is not implemented")

  if(object$order != 1)
    stop("Method for is not implemented for order ", object$order)

  if(length(object$state_vecs) == 0){
    warning("Residuals with State space errors called for model with no time varying effects. Returning NULL")
    return(NULL)
  }

  res = diff(object$state_vecs)

  covs <- vars <- object$state_vars
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

obs_res <- function(object, data, type){
  if(is.null(data) || is.null(object$risk_set))
    stop("Missing risk set or data to compute residuals")

  # Wee need these to check if there is an event in the bin
  new_stop = object$risk_set$stop_new
  new_event = object$risk_set$new_events_flags
  res = list()

  for(i in seq_along(object$risk_set$risk_sets)){
    start_ = object$times[i]
    stop_ = object$times[i + 1]
    r_set = object$risk_set[[1]][[i]]

    tmp_dat = data[r_set, ]
    tmp_dat$tstart = rep(start_, length(r_set))
    tmp_dat$tstop = rep(stop_, length(r_set))

    suppressMessages(p_est <- predict(object, tmp_dat, type = "response", tstart = "tstart", tstop = "tstop")$fits)

    Y = object$risk_set$is_event_in[r_set] == (i - 1)

    if(type == "raw"){
      res[[i]] = cbind(residuals = Y - p_est,
                       p_est = p_est,
                       Y = Y,
                       row_num = r_set)
      next
    }

    res[[i]] = cbind(residuals = (Y - p_est)/sqrt(p_est * (1 - p_est)),
                     p_est = p_est,
                     Y = Y,
                     row_num = r_set)
  }

  return(structure(list(
    residuals = res, type = type), "class" = "fahrmeier_94_res"))
}
