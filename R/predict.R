#' Predict function for the result of \code{\link{ddhazard}}
#' @param object Result of a \code{\link{ddhazard}} call
#' @param new_data New data to base predictions on
#' @param type Either \code{"response"} for predicted probability of death or \code{"term"} for predicted terms in the linear predictor
#' @param tstart Name of the start time column in \code{new_data}. It must corresponds to tstart used in the \code{\link[survival]{Surv}(tstart, tstop, event)} in the \code{formula} passed to \code{\link{ddhazard}}
#' @param tstop same as \code{tstart} for the stop argument
#' @param use_parallel \code{TRUE} if computation for \code{type = "response"} should be computed in parallel with the \code{\link{mcmapply}}. Notice the limitation in the help page of \code{\link{mcmapply}}.
#' @param sds \code{TRUE} if point wise standard deviation should be computed. Convenient if you use functions like \code{\link[splines]{ns}} and you only want one term per term in the right hand site of the \code{formula} used in \code{\link{ddhazard}}
#' @param max_threads Maximum number of threads to use. -1 if it should be determined by a call to \code{\link[parallel]{detectCores}}
#' @param ... Not used
#'
#' @section Term:
#' The result of \code{type = "term"} is a list with the following elements
#' \describe{
#' \item{\code{terms} }{ Is a 3D array. The first dimension is the number of bins, the second dimension is rows in \code{new_data} and the last dimension is the state space terms}
#' \item{\code{sds} }{ Similar to \code{terms} for the point wise confidence intervals using the smoothed co-variance matrices}
#' \item{\code{fixed_terms} }{ Vector of the fixed effect terms for each observation}
#'}
#'
#' @section Response:
#' The result of \code{type = "response"} is a list with the elements below. The function check if there are columns in \code{new_data} which's names match \code{tstart} and \code{tstop}. If not, then each row in new data will get a predicted probability of dying in every bin.
#' \describe{
#' \item{\code{fits} }{ Fitted probability of dying }
#' \item{\code{istart} }{ Vector with the index of the first bin the elements in \code{fits} is in }
#' \item{\code{istop} }{ Vector with the index of the last bin the elements in \code{fits} is in}
#'}
#' @importFrom parallel mcmapply detectCores
#' @export
predict.ddhazard = function(object, new_data,
                                type = c("response", "term"),
                                tstart = "start", tstop = "stop",
                                use_parallel = F, sds = F, max_threads = getOption("ddhazard_max_threads"), ...)
{
  if(!object$model %in% c("logit", exp_model_names))
    stop("Functions for model '", object$model, "' is not implemented")

  type = type[1]
  # TODO: change the code below to use terms object to get the index of
  #       potential fixed effects
  tmp = get_design_matrix(
    data = new_data, response = F, Terms = object$terms, xlev = object$xlev,
    has_fixed_intercept = object$has_fixed_intercept)
  object$formula <- tmp$formula_used

  if(type %in% c("term"))
    return(predict_terms(object, new_data, tmp$X, sds, tmp$fixed_terms))

  if(type %in% c("response"))
    return(predict_response(object, new_data, tmp$X, tstart, tstop, use_parallel, sds, tmp$fixed_terms,
                            max_threads = max_threads))

  stop("Type '", type, "' not implemented in predict.ddhazard")
}

predict_terms <- function(object, new_data, m, sds, fixed_terms){
  # Find the string index maps
  # We have to format the string to a regexp
  term_names_org = c("(Intercept)", attr(object$formula,"term.labels"))
  term_names = stringr::str_replace_all(term_names_org, "(\\W)", "\\\\\\1") # add escape version of charecters

  var_names <- colnames(object$state_vecs)
  if(object$order == 2)
    var_names <- var_names[1:(length(var_names) / 2)]
  terms_to_vars = sapply(term_names, function(t_name) which(grepl(t_name, var_names)))
  found_match <- which(lapply(terms_to_vars, length)  > 0)

  terms_to_vars <- terms_to_vars[found_match]

  stopifnot(!duplicated(unlist(terms_to_vars)))
  stopifnot(length(setdiff(unlist(terms_to_vars), seq_along(var_names))) == 0)

  term_names_org <- term_names_org[found_match]

  # Predict terms
  d <- length(object$times)
  terms_res = array(NA_real_, dim = c(d, nrow(new_data), length(term_names_org)), dimnames = list(NULL, NULL, term_names_org))

  sds_res = if(sds) terms_res else NULL

  for(i in seq_along(term_names_org)){
    terms_res[, , i] = object$state_vecs[ , terms_to_vars[[i]], drop = F] %*%
      t(m[, terms_to_vars[[i]], drop = F])
    if(!sds)
      next

    for(j in seq_len(d))
      sds_res[j, , i] = sqrt(diag(m[, terms_to_vars[[i]], drop = F] %*%
                                    object$state_vars[terms_to_vars[[i]], terms_to_vars[[i]], j] %*%
                                    t(m[, terms_to_vars[[i]], drop = F])))
  }

  fixed_terms <- fixed_terms %*% object$fixed_effects

  return(list(terms = terms_res, sds = sds_res, fixed_terms = fixed_terms))
}

predict_response <- function(
  object, new_data, m, tstart, tstop, use_parallel, sds, fixed_terms,
  max_threads){
  # Change drop behavior inside this function
  old <- `[`
  `[` <- function(...) { old(..., drop=FALSE) }

  # Check order of random walk
  if(object$order > 1)
    warning("Predict not test with new data for order ", object$order)

  # Check if start and stop is provided. If so we need to use these
  # if not, we predict for the sample period
  d <- length(object$times) - 1
  if(all(c(tstart, tstop) %in% colnames(new_data))){
    message("start and stop times ('", tstart, "' and '", tstop, "') are in data. Prediction will match these periods")

    # Find min start. Throw error if before time zero
    start = new_data[[tstart]]
    stop_ = new_data[[tstop]]

    start_order = order(start, method = "radix") - 1L
    start <- round_if_almost_eq(start, start_order, object$times)
    stop_order = order(stop_, method = "radix") - 1L
    stop_ <- round_if_almost_eq(stop_, stop_order, object$times)

  } else{
    message("start and stop times ('", tstart, "' and '", tstop, "') are not in data columns. Each row in new_data will get a row for every bin")
    n_obs <- nrow(m)
    m <- m[sapply(1:nrow(m), rep.int, times = d), ]
    start <- rep(object$times[-(d + 1)], n_obs)
    stop_ <- rep(object$times[-1], n_obs)
  }

  if(min(start) < object$times[1])
    stop("First start time is before time zero")

  # Make prediction of covariates if last stop is beyond last period
  parems = object$state_vecs
  times = object$times

  max_stop = max(stop_)
  max_times = tail(object$times, 1)

  # Check if we have to predict state variables in the future
  if(max_stop > (last_time <- tail(object$times, 1))){
    last_gab = diff(tail(object$times, 2))
    new_times = seq(last_time, max_stop, last_gab)[-1]

    n <- length(new_times)
    if(!isTRUE(all.equal(new_times[n], max_stop))){
      new_times <- c(new_times, new_times[n] + last_gab)
    } else if(new_times[n] < max_stop)
      new_times[n] <- max_stop # needed when we use findInterval later

    n_cols = dim(parems)[2]
    parems = rbind(parems, matrix(NA_real_, nrow = length(new_times), ncol = n_cols))
    if(object$order > 1)
      warning("Currently forecasting wihtout drift from higher than first order effects")
    for(t in seq_along(new_times) + length(times))
      parems[t, ] = parems[t - 1, ]
    times = c(times, new_times)
  }

  if(length(parems) > 0){
    parems = parems[-1, ] # remove first row it is the initial state space vector
    parems = parems[, 1:(dim(parems)[2] / object$order)] # We only need the current estimates (relevant for higher than 1. order)
  }

  # Round if needed. Throw warning if we do so
  int_start = findInterval(start, times)
  if(any(start - times[int_start] > 0) && !object$model %in% exp_model_names)
    warning("Some start times are rounded down")

  int_stop_ = findInterval(stop_, times, left.open = TRUE) # TODO: better way to deal with equality in the stop time?
  if(any(times[int_stop_] - stop_ > 0) && !object$model %in% exp_model_names)
    warning("Some stop times are rounded up")

  # Make function to predict for each observations
  # The function was slow on smaller data sets do to compilation. Thus, the
  # solution below which should only yeild to calls to compile
  .env <- new.env(parent = asNamespace("dynamichazard"))
  .env$discrete_hazard_func = object$discrete_hazard_func
  .env$times <- times
  .env$parems <- parems

  apply_func = with(
    .env, {

      FUN <- function(t, i, i_max, tstart, tstop, x, offset){
        tart <- if(i == 0) tstart else times[t]
        ttop <- if(i == i_max) tstop else times[t + 1]

        discrete_hazard_func(parems[t, ] %*% x + offset, ttop - tart)
      }

      function(istart, istop, x, tstart, tstop, offset){
        #####
        # Compute
        ts = istart:istop
        is = seq_along(ts) - 1
        i_max <- istop - istart
        survival_probs = 1 - mapply(
          FUN,
          t = ts, i = is, i_max = i_max, x = x, tstart =  tstart,
          tstop = tstop, offset = offset)
        1 - prod(survival_probs)
      }})

  # Compute hazard
  args <- list(
    FUN = apply_func,
    x = apply(m, 1, list),
    offset = if(length(object$fixed_effects) == 0)
      rep(0, length(tstart)) else fixed_terms %*% object$fixed_effects,
    istart = int_start, istop = int_stop_, tstart = start, tstop = stop_,
    USE.NAMES = FALSE)

  if(use_parallel){
    no_cores <- detectCores()
    no_cores <- if(is.na(no_cores)){
      1

    } else if (.Platform$OS.type == "windows"){
      1

    } else
      max(min(no_cores - 1, ceiling(nrow(m) / 25)), 1)

    if(max_threads > 0)
      no_cores = min(no_cores, max_threads)

    args <- c(args, list(mc.cores = no_cores))
    fits <- do.call(mcmapply, args)

  }
  else{
    fits <- do.call(mapply, args)

  }

  return(list(
    fits = fits, istart = times[int_start], istop = times[int_stop_]))
}
