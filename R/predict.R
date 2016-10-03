#' Predict function for result of ddhazard
#' @export
predict.fahrmeier_94 = function(object, new_data,
                                type = c("response", "term"),
                                tstart = "start", tstop = "stop",
                                use_parallel = F, sds = F, ...)
{
  if(!object$model %in% c("logit", "exponential"))
    stop("Functions for model '", object$model, "' is not implemented")

  type = type[1]
  m = get_design_matrix(formula = object$formula, data = new_data, response = F)$X

  if(type %in% c("term"))
    return(predict_terms(object, new_data, m, sds))

  if(type %in% c("response"))
    return(predict_response(object, new_data, m, tstart, tstop, use_parallel, sds))

  stop("Type '", type, "' not implemented in predict.fahrmeier_94")
}

predict_terms <- function(object, new_data, m, sds){
  # Find the string index maps
  # We have to format the string to a regexp
  term_names_org = c("(Intercept)", attr(object$formula,"term.labels"))
  term_names = stringr::str_replace_all(term_names_org, "(\\W)", "\\\\\\1") # add escape version of charecters

  var_names = colnames(object$a_t_d_s)
  terms_to_vars = sapply(term_names, function(t_name) which(grepl(t_name, var_names)))

  stopifnot(!duplicated(unlist(terms_to_vars)))
  stopifnot(length(setdiff(unlist(terms_to_vars), seq_along(var_names))) == 0)

  # Predict terms
  d <- length(object$times)
  terms_res = array(NA_real_, dim = c(d, nrow(new_data), length(term_names)), dimnames = list(NULL, NULL, term_names_org))

  sds_res = if(sds) terms_res else NULL

  for(i in seq_along(term_names)){
    terms_res[, , i] = object$a_t_d_s[ , terms_to_vars[[i]], drop = F] %*%
      t(m[, terms_to_vars[[i]], drop = F])
    if(!sds)
      next

    for(j in seq_len(d))
      sds_res[j, , i] = sqrt(diag(m[, terms_to_vars[[i]], drop = F] %*%
                                    object$V_t_d_s[j, terms_to_vars[[i]], terms_to_vars[[i]]] %*%
                                    t(m[, terms_to_vars[[i]], drop = F])))
  }

  return(list(terms = terms_res, sds = sds_res))
}

predict_response <- function(object, new_data, m, tstart, tstop, use_parallel, sds){
  # Check order of random walk
  if(object$order > 1)
    warning("Predict not test with new data for order ", object$order)

  # Check if start and stop is provided. If so we need to use these
  # if not, we predict for the sample period
  d <- length(object$times) - 1
  if(all(c(tstart, tstop) %in% colnames(new_data))){
    message("start and stop times ('", tstart, "' and '", tstop, "') are in data. Prediction will match these periods")

    # Find min start. Throw error if before time zero
    start = new_data[, tstart]
    stop_ = new_data[, tstop]
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
  parems = object$a_t_d_s
  times = object$times

  max_stop = max(stop_)
  max_times = tail(object$times, 1)

  # Check if we have to predict state variables in the future
  if(max_stop > tail(object$times, 1)){
    last_gab = diff(tail(object$times, 2))
    new_times = max_times + last_gab*(1:ceiling((max_stop - max_times) / last_gab))

    n_cols = dim(parems)[2]
    parems = rbind(parems, matrix(NA_real_, nrow = length(new_times), ncol = n_cols))
    if(object$order > 1)
      warning("Currently forecasting wihtout drift from higher than first order effects")
    for(t in seq_along(new_times) + length(times))
      parems[t, ] = parems[t - 1, ]
    times = c(times, new_times)
  }

  parems = parems[-1, ] # remove first row it is the initial state space vector
  parems = parems[, 1:(dim(parems)[2] / object$order)] # We only need the current estimates (relevant for higher than 1. order)

  # Round if needed. Throw error if so
  int_start = findInterval(start, times)
  if(any(start - times[int_start] > 0) && object$model != "exponential")
    warning("Some start times are rounded down")

  int_stop_ = findInterval(stop_, times, left.open = T) # TODO: better way to deal with equality in the stop time?
  if(any(times[int_stop_] - stop_ > 0) && object$model != "exponential")
    warning("Some stop times are rounded up")

  # Make function to predict for each observations
  # assume that covariates do not change
  hazard_func = object$hazard_func
  tmp_func = function(x_, istart, istop, tstart, tstop){

    i <- 0
    i_max <- istop - istart

    survival_probs = 1 - sapply(istart:istop, function(t){
      tart <- if(i == 0) tstart else times[t]
      ttop <- if(i == i_max) tstop else times[t + 1]
      i <<- i + 1
      hazard_func(parems[t, ] %*% x_, tstart = tart, tstop = ttop)
    })

    1 - prod(survival_probs)
  }

  # Compute hazard
  if(use_parallel){
    tryCatch({
      no_cores <- min(parallel::detectCores() - 1, ceiling(nrow(m) / 25))
      cl <- parallel::makeCluster(no_cores)
      parallel::clusterExport(cl, c("parems", "hazard_func", "times"),
                              envir = environment())

      fits = parallel::parRapply(cl = cl, data.frame(istart = int_start, istop = int_stop_,
                                                     tstart = start, tstop = stop_, x_ = m),
                                 function(row_){
                                   tmp_func(x_ = row_[-(1:4)], istart = row_[1], istop = row_[2],
                                            tstart =  row_[3], tstop =  row_[4])
                                 })

    }, finally = { parallel::stopCluster(cl)})
  }
  else{
    fits = apply(data.frame(istart = int_start, istop = int_stop_,
                            tstart = start, tstop = stop_, x_ = m), 1,
                 function(row_){
                   tmp_func(x_ = row_[-(1:4)], istart = row_[1], istop = row_[2],
                            tstart =  row_[3], tstop =  row_[4])
                 })
  }

  return(list(fits = fits, istart = times[int_start], istop = times[int_stop_]))
}
