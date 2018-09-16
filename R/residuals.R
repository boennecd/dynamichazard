#' @title Residuals Method for ddhazard Object
#' @description Residuals method for the result of a \code{\link{ddhazard}} call.
#'
#' @param object result of \code{\link{ddhazard}} call.
#' @param type type of residuals. Four possible values: \code{"std_space_error"}, \code{"space_error"}, \code{"pearson"} and \code{"raw"}. See the sections below for details.
#' @param data \code{data.frame} with data for the Pearson or raw residuals. This is only needed if the data set is not saved with the \code{object}. Must be the same data set used in the initial call to \code{\link{ddhazard}}.
#' @param ... not used.
#'
#' @section Pearson and raw residuals:
#' Is the result of a call with a \code{type} argument of either \code{"pearson"} or \code{"raw"} for Pearson residuals or raw residuals. Returns a list with class \code{"ddhazard_residual"} with the following elements.
#' \describe{
#' \item{\code{residuals}}{list of residuals for each bin. Each element of the list contains a 2D array where the rows corresponds to the passed \code{data} and columns are the residuals (\code{residuals}), estimated probability of death (\code{p_est}), outcome (\code{Y}) and row number in the initial data set (\code{row_num}). The \code{data} rows will only have a residuals in a given risk list if they are at risk in that risk set.}
#' \item{\code{type}}{the type of residual.}
#'}
#'
#' @section State space errors:
#' Is the result of a call with a \code{type} argument of either \code{"std_space_error"} or \code{"space_error"}. The former is for standardized residuals while the latter is non-standardized. Returns a list with class. \code{"ddhazard_space_errors"} with the following elements:
#' \describe{
#' \item{\code{residuals}}{2D array with either standardized or non-standardized state space errors. The row are bins and the columns are the parameters in the regression.}
#' \item{\code{standardize}}{\code{TRUE} if standardized state space errors.}
#' \item{\code{Covariances}}{3D array with the smoothed co-variance matrix for each set of the state space errors.}
#'}
#'
#' @examples
#'library(dynamichazard)
#'fit <- ddhazard(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
#'  Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
#'  control = ddhazard_control(method = "GMA"))
#'resids <- residuals(fit, type = "pearson")$residuals
#'head(resids[[1]])
#'head(resids[[2]])
#'
#' @export
residuals.ddhazard = function(
  object, type = c("std_space_error", "space_error", "pearson", "raw"),
  data = NULL, ...){
  type = type[1]

  if(!object$model %in% c("logit", exp_model_names))
    stop("Functions for model '",  object$model, "' is not implemented")

  if(type %in% c("std_space_error", "space_error"))
    return(space_errors(object, data, type == "std_space_error"))

  if(type == "pearson" || type == "raw"){
    if(type == "pearson" & !object$model %in% c("logit"))
      stop("Pearsons residuals is not implemented for model '", object$model, "'")

    return(
      obs_res(object, if(is.null(object$data)) data else object$data, type))
  }

  stop("Method '", type, "' not implemented for residuals method")
}

space_errors <- function(object, data, standardize){
  if(!object$method %in% c("EKF"))
    stop("Functions for with method '", object$method, "' is not implemented")

  if(!object$model %in% c("logit"))
    stop("Functions for with model '", object$model, "' is not implemented")

  if(object$order != 1)
    stop("Method for is not implemented for order ", object$order)

  if(length(object$state_vecs) == 0){
    warning("Residuals with State space errors called for model with no time varying effects. Returning NULL")
    return(NULL)
  }

  res = diff(object$state_vecs)

  covs <- vars <- object$state_vars
  for(i in seq_len(dim(vars)[3] - 1))
    covs[, , i] = vars[, , i + 1] + vars[, , i] -
      object$lag_one_cov[, , i] - t(object$lag_one_cov[, , i])

  if(standardize)
    for(i in seq_len(dim(vars)[3] - 1))
      res[i, ] = res[i, ] %*% solve(chol(covs[, , i]))

  return(structure(list(residuals = res,
                        standardize = standardize,
                        Covariances = covs),
                   "class" = "ddhazard_space_errors"))
}

obs_res <- function(object, data, type){
  if(is.null(data) || is.null(object$risk_set))
    stop("Missing risk set or data to compute residuals")

  # we need these to check if there is an event in the bin
  new_stop = object$risk_set$stop_new
  new_event = object$risk_set$new_events_flags
  res = list()

  # get the reponse object
  response <- model.frame(update(object$formula, .~-1), data)[[1]]
  if(attr(response, "type") == "right")
    response <- cbind(rep(0, nrow(response)), response)

  # find the rows tstart and tstop
  tstart <- response[, 1]
  tstop <- response[, 2]

  for(i in seq_along(object$risk_set$risk_sets)){
    start_ = object$times[i]
    stop_ = object$times[i + 1]
    r_set = object$risk_set[[1]][[i]]

    tmp_dat = data[r_set, ]
    tmp_dat$tstart = pmax(tstart[r_set], start_)
    tmp_dat$tstop = pmin(tstop[r_set], stop_)

    suppressMessages(p_est <- predict(
      object, tmp_dat, type = "response", tstart = "tstart",
      tstop = "tstop")$fits)

    Y = object$risk_set$is_event_in[r_set] == (i - 1)

    if(type == "raw"){
      res[[i]] = cbind(
        residuals = Y - p_est, p_est = p_est, Y = Y, row_num = r_set)
      next
    }

    res[[i]] = cbind(
      residuals = (Y - p_est)/sqrt(p_est * (1 - p_est)), p_est = p_est,
      Y = Y, row_num = r_set)
  }

  return(structure(list(
    residuals = res, type = type), "class" = "ddhazard_residual"))
}
