#' Predict Method for ddhazard Object
#'
#' @description Predict method for \code{\link{ddhazard}}.
#'
#' @param object result of a \code{\link{ddhazard}} call.
#' @param new_data new data to base predictions on.
#' @param type either \code{"response"} for predicted probability of an event or \code{"term"} for predicted terms in the linear predictor.
#' @param tstart name of the start time column in \code{new_data}. It must be on the same time scale as the \code{tstart} used in the \code{\link[survival]{Surv}(tstart, tstop, event)} in the \code{formula} passed to \code{\link{ddhazard}}.
#' @param tstop same as \code{tstart} for the stop argument.
#' @param use_parallel not longer supported.
#' @param sds \code{TRUE} if point wise standard deviation should be computed. Convenient if you use functions like \code{\link[splines]{ns}} and you only want one term per term in the right hand site of the \code{formula} used in \code{\link{ddhazard}}.
#' @param max_threads not longer supported.
#' @param ... not used.
#'
#' @details
#' The function check if there are columns in \code{new_data} which names match
#' \code{tstart} and \code{tstop}. If matched, then the bins are found which
#' the start time to the stop time are in. If \code{tstart} and \code{tstop} are not
#' matched then all the bins used in the estimation method will be used.
#'
#' @section Term:
#' The result with \code{type = "term"} is a lists of list each having length
#' equal to \code{nrow(new_data)}. The lists are
#' \describe{
#' \item{\code{terms}}{It's elements are matrices where the first dimension is
#' time and the second dimension is the terms.}
#' \item{\code{sds}}{similar to \code{terms} for the point-wise confidence
#' intervals using the smoothed co-variance matrices. Only added if
#' \code{sds = TRUE}.}
#' \item{\code{fixed_terms}}{contains the fixed (non-time-varying) effect.}
#' \item{\code{varcov}}{similar to \code{sds} but differs by containing the whole
#' covariance matrix for the terms. It is a 3D array where the third dimension is
#' time. Only added if \code{sds = TRUE}.}
#' \item{\code{start}}{numeric vector with start time for each time-varying term.}
#' \item{\code{tstop}}{numeric vector with stop time for each time-varying term.}
#'}
#'
#' @section Response:
#' The result with \code{type = "response"} is a list with the elements below.
#' If \code{tstart} and \code{tstop} are matched in columns in \code{new_data},
#' then the probability will be for having an event between \code{tstart} and \code{tstop}
#' conditional on no events before \code{tstart}.
#' \describe{
#' \item{\code{fits}}{fitted probability of an event.}
#' \item{\code{istart}}{numeric vector with start time for each element in \code{fits}.}
#' \item{\code{istop}}{numeric vector with stop time for each element in \code{fits}.}
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
#' # probability of an event between time 0 and 2000 with bili = 3
#' predict(fit, type = "response", new_data =
#'           data.frame(time = 0, status = 2, bili = 3, tstart = 0, tstop = 2000),
#'         tstart = "tstart", tstop = "tstop")
#'
#' @importFrom parallel mcmapply detectCores
#' @export
predict.ddhazard = function(object, new_data,
                                type = c("response", "term"),
                                tstart = "start", tstop = "stop",
                                use_parallel, sds = F, max_threads, ...)
{
  if(!object$model %in% c("logit", "cloglog", exp_model_names))
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
      if(object$order > 1)
        stop(sQuote("sds"), " = TRUE is not implemented with ",
             sQuote("order"), " > 1")

      new_state_vars <- array(dim = c(dim(state_vars)[1:2], nrow(params)))
      new_state_vars[, , 1:dim(state_vars)[3]] <- state_vars

      Q_t <- Q. * last_gab
      for(i in seq_along(new_times) + length(otimes))
        new_state_vars[, , i] <-
          tcrossprod(F_ %*% new_state_vars[, , i - 1L], F_) + Q_t

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

  Jt <- matrix(0., ncol(m), length(term_names_org))
  for(i in seq_along(term_names_org))
    for(j in terms_to_vars[i])
      Jt[j, i] <- 1
  XJ <- if(length(m) == 0)
    replicate(nrow(m), list(m)) else lapply(split(m, row(m)), "*", y = Jt)

  out <- mapply(
    .predict_terms, sta = start, sto = stop_, XJ = XJ,
    ista = int_start, isto = int_stop_, ft = fixed_terms,
    MoreArgs = list(
      params = params, state_vars = state_vars, comp_sds = sds, byl =
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

.predict_terms <- function(XJ, sta, sto, ista, isto, ft, comp_sds, otimes,
                           params, state_vars, byl)
{
  ts = ista:isto
  is = seq_along(ts) - 1L
  i_max <- isto - ista

  # compute terms
  O <- if(length(XJ) == 0)
    params[ts, , drop = FALSE] else params[ts, ] %*% XJ

  # compute covariance matrix
  if(comp_sds)
  {
    q <- ncol(XJ)
    d <- length(ts)
    varcov <- apply(state_vars[, , ts, drop = FALSE], 3, function(x)
      crossprod(XJ, x %*% XJ))
    dim(varcov) <- c(q, q, d)

    sds_res <- sqrt(apply(varcov, 3, diag))
    sds_res <- if(is.vector(sds_res))
      as.matrix(sds_res) else t(sds_res)

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

predict_response <- function(
  object, new_data, m, tstart, tstop, fixed_terms, keep_separate = FALSE){
  trs <- predict_terms(
    object = object, new_data = new_data, m = m, fixed_terms = fixed_terms,
    tstart = tstart, tstop = tstop, sds = FALSE)

  has_times <- !keep_separate && all(c(tstart, tstop) %in% colnames(new_data))
  pevent <- with(
    trs, mapply(
      .predict_response, terms. = terms, fixed_terms = fixed_terms,
      tstart = tstart, tstop = tstop, MoreArgs = list(
        object = object, has_times = has_times, ddfunc =
          object$discrete_hazard_func)))
  if(has_times){
    return(list(
      fits = pevent, istart = sapply(trs$tstart, "[[", 1),
      istop = sapply(trs$tstop, function(x) x[length(x)])))
  }

  list(fits = drop(pevent), istart = unlist(trs$tstart),
       istop = unlist(trs$tstop))
}

.predict_response <- function(
  terms., fixed_terms, tstart, tstop, object, has_times, ddfunc)
  {
    probs <- ddfunc(
      eta = rowSums(terms.) + fixed_terms, at_risk_length = tstop - tstart)

    if(has_times)
      return(1  - prod(1 - probs))

    probs
}

#' @name ddsurvcurve
#' @title Create and plot survival curves
#' @description
#' The function creates a predicted survival curve for a new observation using
#' a estimated \code{ddhazard} model from \code{\link{ddhazard}}. The predicted
#' curve is based on the predicted mean path of the state vector. Thus, the
#' survival curve will not be a "mean" curve due to the non-linear relation between
#' the probability of an event and the state vector.
#'
#' @inheritParams predict.ddhazard
#' @param object a \code{ddhazard} object.
#' @param new_data a \code{data.frame} with the new data for the observation
#' who the survival curve should be for. It can have more rows if \code{tstart}
#' and \code{tstop} is supplied. The rows need to be consecutive and non-overlapping
#' time intervals.
#' @param x a \code{ddsurvcurve} object.
#' @param y not used.
#' @param xlab \code{xlab} passed to \code{plot}.
#' @param ylab \code{ylab} passed to \code{plot}.
#' @param ylim \code{ylim} passed to \code{plot}.
#' @param xaxs \code{xaxs} passed to \code{plot}.
#' @param yaxs \code{yaxs} passed to \code{plot}.
#' @param col \code{col} passed to \code{lines}.
#' @param lty \code{lty} passed to \code{lines}.
#' @param lwd \code{lwd} passed to \code{lines}.
#'
#' @return
#' \code{ddsurvcurve} returns an object of class \code{ddsurvcurve}. It elements are the predicted
#' discrete  survival curve, time points for the survival curve, point of the first
#' time period, the call, the discrete probabilities of an event in each interval
#' conditional on survival up to that point, and the name of the distribution
#' family. It should be seen as a plug-in estimate.
#'
#' @seealso
#' \code{\link{ddhazard}}, and \code{\link{predict.ddhazard}}.
#'
#' @examples
#' #####
#' # example with continuous time model
#' # setup data set. See `vignette("timedep", "survival")`
#' library(dynamichazard)
#' temp <- subset(pbc, id <= 312, select=c(id:sex, stage))
#' pbc2 <- tmerge(temp, temp, id=id, death = event(time, status))
#' pbc2 <- tmerge(pbc2, pbcseq, id = id, bili = tdc(day, bili))
#'
#' # fit model
#' f1 <- ddhazard(
#'   Surv(tstart, tstop, death == 2) ~ ddFixed(log(bili)), pbc2, id = pbc2$id,
#'   max_T = 3600, Q_0 = 1, Q = 1e-2, by = 100, model = "exponential",
#'   control = ddhazard_control(method = "EKF", eps = 1e-4, n_max = 1000,
#'                              fixed_terms_method = "M_step"))
#'
#' # predict with default which is all covariates set to zero
#' ddcurve <- ddsurvcurve(f1)
#' par(mar = c(4.5, 4, .5, .5))
#' plot(ddcurve, col = "DarkBlue", lwd = 2)
#'
#' # compare with cox model
#' f2 <- coxph(Surv(tstart, tstop, death == 2) ~ log(bili), data = pbc2)
#' nw <- data.frame(bili = 1, tstart = 0, tstop = 3000)
#' lines(survfit(f2, newdata = nw))
#'
#' # same as above but with bili = 3
#' nw <- data.frame(bili = 3)
#' lines(ddsurvcurve(f1, new_data = nw), col = "DarkBlue")
#' lines(survfit(f2, newdata = nw))
#'
#' # change to time-varying slope
#' f3 <- ddhazard(
#'   Surv(tstart, tstop, death == 2) ~ log(bili), pbc2, id = pbc2$id,
#'   max_T = 3600, Q_0 = diag(1, 2), Q = diag(1e-2, 2), by = 100, model = "exponential",
#'   control = ddhazard_control(method = "EKF", eps = 1e-4, n_max = 1000))
#'
#' # example with time-varying coefficient
#' nw <- data.frame(
#'   bili = c(2.1, 1.9, 3.3, 3.9, 3.8, 3.6, 4, 4.9, 4.2, 5.7, 10.2),
#'   tstart = c(0L, 225L, 407L, 750L, 1122L, 1479L, 1849L, 2193L, 2564L, 2913L,
#'              3284L),
#'   tstop = c(225L, 407L, 750L, 1122L, 1479L, 1849L, 2193L, 2564L, 2913L,
#'             3284L, 3600L))
#' ddcurve <- ddsurvcurve(f3, new_data = nw, tstart = "tstart", tstop = "tstop")
#' lines(ddcurve, "darkorange", lwd = 2)
#'
#' # can condition on survival up to some time
#' ddcurve <- ddsurvcurve(f3, new_data = nw[-(1:5), ], tstart = "tstart",
#'                        tstop = "tstop")
#' lines(ddcurve, lty = 2, lwd = 2)
#'
#' #####
#' # example with discrete time model
#' # head-and-neck cancer study data. See Efron, B. (1988) doi:10.2307/2288857
#' is_censored <- c(
#'   6, 27, 34, 36, 42, 46, 48:51, 51 + c(15, 30:28, 33, 35:37, 39, 40, 42:45))
#' head_neck_cancer <- data.frame(
#'   id = 1:96,
#'   stop = c(
#'     1, 2, 2, rep(3, 6), 4, 4, rep(5, 8),
#'     rep(6, 7), 7, 8, 8, 8, 9, 9, 10, 10, 10, 11, 14, 14, 14, 15, 18, 18, 20,
#'     20, 37, 37, 38, 41, 45, 47, 47,
#'     2, 2, 3, rep(4, 4), rep(5, 5), rep(6, 5),
#'     7, 7, 7, 9, 10, 11, 12, 15, 16, 18, 18, 18, 21,
#'     21, 24, 25, 27, 36, 41, 44, 52, 54, 59, 59, 63, 67, 71, 76),
#'   event = !(1:96 %in% is_censored),
#'   group = factor(c(rep(1, 45 + 6), rep(2, 45))))
#'
#' # fit model
#' h1 <- ddhazard(
#'   Surv(stop, event) ~ group, head_neck_cancer, by = 1, max_T = 45,
#'   Q_0 = diag(2^2, 2), Q = diag(.01^2, 2), control = ddhazard_control(
#'     method = "GMA", eps = 1e-4, n_max = 200))
#'
#' # plot predicted survival curve. Notice the steps since the model is discrete
#' nw <- data.frame(group = factor(1, levels = 1:2), tstart = 0, tstop = 30)
#' ddcurve <- ddsurvcurve(h1, new_data = nw, tstart = "tstart",
#'                        tstop = "tstop")
#' plot(ddcurve, col = "Darkblue")
#'
#' nw$group <- factor(2, levels = 1:2)
#' ddcurve <- ddsurvcurve(h1, new_data = nw, tstart = "tstart",
#'                        tstop = "tstop")
#' lines(ddcurve, col = "DarkOrange")
#'
#' # compare with KM
#' lines(survfit(Surv(stop, event) ~ 1, head_neck_cancer, subset = group == 1),
#'       col = "DarkBlue")
#' lines(survfit(Surv(stop, event) ~ 1, head_neck_cancer, subset = group == 2),
#'       col = "DarkOrange")
#'
#' @export
#' @importFrom utils head
ddsurvcurve <- function(object, new_data, tstart = "", tstop = ""){
  #####
  # find predicted probabilities
  if(missing(new_data)){
    if(!"(Intercept)" %in% colnames(object$state_vecs))
      stop(sQuote("ddsurvcurve"), " not implemented for model without intercept",
           " and call without ", sQuote("new_data"))

    ti <- object$times
    trs <- object$discrete_hazard_func(
      object$state_vecs[-1, "(Intercept)"], at_risk_length = diff(object$times))
    trs <- list(fits = trs, istart = head(ti, -1), istop = ti[-1])

  } else {
    # TODO: change the code below to use terms object which has been added in
    #       in version 0.5.0 to get the index of potential fixed effects
    tmp = get_design_matrix(
      formula = object$formula,
      data = new_data, response = F, Terms = object$terms, xlev = object$xlev,
      has_fixed_intercept = object$has_fixed_intercept)
    object$formula <- tmp$formula_used
    trs <- predict_response(
      object = object, new_data = new_data, m = tmp$X,
      tstart = tstart, tstop = tstop, fixed_terms = tmp$fixed_terms,
      keep_separate = TRUE)

  }

  # check if there are no non-consecutive or overlapping intervals
  psurv <- cumprod(1 - unlist(trs$fits))
  sta <- unlist(trs$istart)
  sto <- unlist(trs$istop)
  if(!isTRUE(all.equal(sta[-1], head(sto, -1))))
    stop("non-consecutive periods is passed")

  structure(
    list(psurv = psurv, time = sto, start = sta[1], call = match.call(),
         dhazard = unlist(trs$fits), family = object$family$name()),
    class = "ddsurvcurve")
}
