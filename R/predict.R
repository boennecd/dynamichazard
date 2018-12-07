#' Predict Method for ddhazard Object
#'
#' @description Predict method for \code{\link{ddhazard}}.
#'
#' @param object result of a \code{\link{ddhazard}} call.
#' @param new_data new data to base predictions on.
#' @param type either \code{"response"} for predicted probability of death or \code{"term"} for predicted terms in the linear predictor.
#' @param tstart name of the start time column in \code{new_data}. It must be on the same time scale as the tstart used in the \code{\link[survival]{Surv}(tstart, tstop, event)} in the \code{formula} passed to \code{\link{ddhazard}}.
#' @param tstop same as \code{tstart} for the stop argument.
#' @param use_parallel \code{TRUE} if computation for \code{type = "response"} should be computed in parallel with the \code{\link{mcmapply}}. Notice the limitation in the help page of \code{\link{mcmapply}}.
#' @param sds \code{TRUE} if point wise standard deviation should be computed. Convenient if you use functions like \code{\link[splines]{ns}} and you only want one term per term in the right hand site of the \code{formula} used in \code{\link{ddhazard}}.
#' @param max_threads maximum number of threads to use. -1 if it should be determined by a call to \code{\link[parallel]{detectCores}}.
#' @param ... not used.
#'
#' @section Term:
#' The result with \code{type = "term"} is a list with the following elements
#' \describe{
#' \item{\code{terms}}{is a 3D array. The first dimension is the number of bins, the second dimension is rows in \code{new_data} and the last dimension is the state space terms.}
#' \item{\code{sds}}{similar to \code{terms} for the point wise confidence intervals using the smoothed co-variance matrices.}
#' \item{\code{fixed_terms}}{vector of the fixed effect terms for each observation.}
#'}
#'
#' @section Response:
#' The result with \code{type = "response"} is a list with the elements below. The function check if there are columns in \code{new_data} which names match \code{tstart} and \code{tstop}. If not, then each row in new data will get a predicted probability of dying in every bin.
#' \describe{
#' \item{\code{fits}}{fitted probability of dying.}
#' \item{\code{istart}}{vector with the index of the first bin the elements in \code{fits} is in.}
#' \item{\code{istop}}{vector with the index of the last bin the elements in \code{fits} is in.}
#'}
#'
#' @examples
#' fit <- ddhazard(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
#'  Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
#'  control = ddhazard_control(method = "GMA"))
#' predict(fit, type = "response", new_data =
#'  data.frame(time = 0, status = 2, bili = 3))
#' predict(fit, type = "term", new_data =
#'  data.frame(time = 0, status = 2, bili = 3))
#'
#' @importFrom parallel mcmapply detectCores
#' @export
predict.ddhazard = function(object, new_data,
                                type = c("response", "term"),
                                tstart = "start", tstop = "stop",
                                use_parallel, sds = F, max_threads, ...)
{
  if(!object$model %in% c("logit", exp_model_names))
    stop("Functions for model '", object$model, "' is not implemented")
  if(!missing(max_threads) || !missing(use_parallel))
    warning(sQuote("max_threads"), " and ", sQuote("use_parallel"),
            " is not used anymore")

  type = type[1]
  # TODO: change the code below to use terms object which has been added in
  #       in version 0.5.0 to get the index of potential fixed effects
  tmp = get_design_matrix(
    formula = object$formula,
    data = new_data, response = F, Terms = object$terms, xlev = object$xlev,
    has_fixed_intercept = object$has_fixed_intercept)
  object$formula <- tmp$formula_used

  if(type %in% c("term"))
    return(predict_terms(
      object, new_data, tmp$X, tmp$fixed_terms, tstart, tstop, sds = sds))

  if(type %in% c("response"))
    return(predict_response(
      object = object, new_data = new_data, m = tmp$X,
      tstart = tstart, tstop = tstop, fixed_terms = tmp$fixed_terms))

  # object, new_data, m, tstart, tstop, fixed_terms

  stop("Type '", type, "' not implemented in predict.ddhazard")
}

predict_terms <- function(object, new_data, m, sds, fixed_terms, tstart,
                          tstop){
  #####
  # find map to terms
  # we have to format the string to a regexp due as we do not have the
  # `model.matrix` output...
  term_names_org = c("(Intercept)", attr(object$formula, "term.labels"))
  # add escape version of charecters
  term_names = gsub("(\\W)", "\\\\\\1", term_names_org, perl = TRUE)

  var_names <- colnames(object$state_vecs)
  if(object$order == 2)
    var_names <- var_names[1:(length(var_names) / 2)]
  terms_to_vars = sapply(
    term_names, function(t_name) which(grepl(t_name, var_names)))
  found_match <- which(lapply(terms_to_vars, length)  > 0)

  terms_to_vars <- terms_to_vars[found_match]

  stopifnot(!duplicated(unlist(terms_to_vars)))
  stopifnot(length(setdiff(unlist(terms_to_vars), seq_along(var_names))) == 0)

  term_names_org <- term_names_org[found_match]

  #####
  # find timer periods
  otimes <- object$times
  d <- length(otimes) - 1L
  if(all(c(tstart, tstop) %in% colnames(new_data))){
    message("start and stop times ('", tstart, "' and '", tstop,
            "') are in data. Prediction will match these periods")

    # Find min start. Throw error if before time zero
    start = new_data[[tstart]]
    stop_ = new_data[[tstop]]

    start_order = order(start, method = "radix") - 1L
    start <- round_if_almost_eq(start, start_order, otimes)
    stop_order = order(stop_, method = "radix") - 1L
    stop_ <- round_if_almost_eq(stop_, stop_order, otimes)

  } else{
    message(
      "start and stop times ('", tstart, "' and '", tstop,
      "') are not in data columns. Each row in new_data will get a row for every bin")
    n_obs <- nrow(m)
    start <- rep(otimes[1L]    , n_obs)
    stop_ <- rep(otimes[d + 1L], n_obs)

  }

  if(min(start) < otimes[1])
    stop("First start time is before time zero")

  #####
  # make prediction of covariates if last stop is beyond last period
  params = object$state_vecs
  F_ <- object$F_
  Q. <- object$Q
  state_vars <- object$state_vars

  max_stop = max(stop_)
  max_times = tail(otimes, 1)

  # check if we have to predict state variables in the future
  if(max_stop > (last_time <- tail(otimes, 1))){
    last_gab = diff(tail(otimes, 2))
    new_times = seq(last_time, max_stop, last_gab)[-1]

    n <- length(new_times)
    if(!isTRUE(all.equal(new_times[n], max_stop))){
      new_times <- c(new_times, new_times[n] + last_gab)

    } else if(new_times[n] < max_stop)
      new_times[n] <- max_stop # needed when we use findInterval later

    # predict state and covariance matrix of prediction
    params = rbind(params, matrix(
      NA_real_, nrow = length(new_times), ncol = dim(params)[2]))

    for(i in seq_along(new_times) + length(otimes))
      params[i, ] = F_ %*% params[i - 1, ]

    if(sds){
      if(order > 1)
        stop(sQuote("sds"), " = TRUE is not implemented with ",
             sQuote("order"), " > 1")

      new_state_vars <- array(dim = c(dim(state_vars)[1:2], nrow(params)))
      new_state_vars[, , 1:dim(state_vars)[3]] <- state_vars

      for(i in seq_along(new_times) + length(times))
        new_state_vars[, , i] <-
          tcrossprod(F_ %*% new_state_vars[, , i - 1L], F_) + Q.

      state_vars <- new_state_vars
    }

    otimes <- c(otimes, new_times)
  }

  if(length(params) > 0){
    # remove first row it is the initial state vector
    # we only need the current estimates (relevant for higher than 1. order)
    keep <- 1:(dim(params)[2] / object$order)
    params <- params[-1, keep, drop = FALSE]

    state_vars <- state_vars[keep, keep, -1, drop = FALSE]
  }

  #####
  # round if needed. Post warning if we do so
  int_start = findInterval(start, otimes)
  if(any(start - otimes[int_start] > 0) && !object$model %in% exp_model_names)
    warning("Some start times are rounded down")

  int_stop_ = findInterval(stop_, otimes, left.open = TRUE)
  if(any(otimes[int_stop_] - stop_ > 0) && !object$model %in% exp_model_names)
    warning("Some stop times are rounde")

  #####
  # predict term
  fixed_terms <- fixed_terms %*% object$fixed_effects
  out <- mapply(
    .predict_terms, sta = start, sto = stop_, x = apply(m, 1, list),
    ista = int_start, isto = int_stop_, ft = fixed_terms,
    MoreArgs = list(
      params = params, state_vars = state_vars, comp_sds = sds,
      term_idx = terms_to_vars, term_na = term_names_org, byl =
        diff(tail(otimes, 2)), otimes = otimes), SIMPLIFY = FALSE)

  # predict terms
  list(
    terms = lapply(out, "[[", "terms"),
    sds = if(sds) lapply(out, "[[", "sds") else NULL,
    fixed_terms = lapply(out, "[[", "fixed_terms"),
    varcov = if(sds) lapply(out, "[[", "varcov") else NULL,
    tstart = lapply(out, "[[", "tstart"),
    tstop = lapply(out, "[[", "tstop"))
}

.predict_terms <- function(x, sta, sto, ista, isto, ft, comp_sds, otimes,
                           term_idx, term_na, params, state_vars, byl)
{
  ts = ista:isto
  is = seq_along(ts) - 1L
  i_max <- isto - ista

  # compute terms
  x <- x[[1]]
  J <- matrix(0., length(x), length(term_idx))
  for(i in seq_along(term_na))
    for(j in term_idx[i])
      J[j, i] <- 1
  XJ <- diag(x) %*% J
  O <- params[ts, ] %*% XJ

  # compute covariance matrix
  if(comp_sds)
  {
    varcov <- array(dim = c(length(term_na), length(term_na), length(ts)))
    for(i in 1:length(ts))
      varcov[, , i] <- crossprod(XJ, state_vars[, , ts[i]] %*% XJ)

    sds_res <- sqrt(t(apply(varcov, 3, diag)))

  } else {
    varcov  <- NULL
    sds_res <- NULL
  }

  # start and stop times
  stas <- otimes[ts]
  stas[1] <- sta
  stos <- otimes[ts + 1L]
  stos[length(stos)] <- sto

  # return
  list(terms = O, sds = sds_res, fixed_terms = rep(ft, length(ts)),
       varcov = varcov, tstart = stas, tstop = stos)
}

#' @importFrom utils tail
predict_response <- function(
  object, new_data, m, tstart, tstop, fixed_terms){
  trs <- predict_terms(
    object = object, new_data = new_data, m = m, fixed_terms = fixed_terms,
    tstart = tstart, tstop = tstop, sds = FALSE)

  has_times <- all(c(tstart, tstop) %in% colnames(new_data))
  pevent <- with(
    trs, mapply(
      .predict_response, terms. = terms, fixed_terms = fixed_terms,
      tstart = tstart, tstop = tstop, MoreArgs = list(object = object,
                                                      has_times = has_times)))
  if(has_times){
    return(list(
      fits = pevent, istart = sapply(trs$tstart, "[[", 1),
      istop = sapply(trs$tstop, function(x) x[length(x)])))
  }

  list(fits = drop(pevent), istart = unlist(trs$tstart),
       istop = unlist(trs$tstop))
}

.predict_response <- function(
  terms., fixed_terms, tstart, tstop, object, has_times)
  {
    probs <- object$discrete_hazard_func(
      eta = rowSums(terms.) + fixed_terms, at_risk_length = tstop - tstart)

    if(has_times)
      return(1  - prod(1 - probs))

    probs
  }
