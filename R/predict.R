#' Predict function for result of ddhazard
#' @export
predict.fahrmeier_94 = function(object, new_data, level = pnorm(1.96/2), type = "not_term",
                                method = "delta", tstart = "start", tstop = "stop",
                                use_parallel = F, sds = F, ...)
{
  m = benssurvutils::get_design_matrix(formula = object$formula, data = new_data, response = F)$X

  if(type %in% c("term", "term_var_diff")){
    if(length(object$risk_set) <= 1 && is.na(object$risk_set))
      stop("Terms requires the risk set")

    # Find the string index maps
    # We have to format the string to a regexp
    term_names_org = c("(Intercept)", attr(object$formula,"term.labels"))
    term_names = gsub("\\[", "\\\\[", term_names_org, perl = T)
    term_names = gsub("\\]", "\\\\]", term_names, perl = T)
    term_names = gsub('\\"', '\\\\"', term_names, perl = T)
    term_names = gsub('\\(', '\\\\(', term_names, perl = T)
    term_names = gsub('\\)', '\\\\)', term_names, perl = T)
    term_names = paste("^", term_names, sep = "")

    var_names = colnames(object$a_t_d_s)

    terms_to_vars = sapply(term_names, function(t_name) which(grepl(t_name, var_names)))
    stopifnot(!duplicated(unlist(terms_to_vars)))
    stopifnot(length(setdiff(unlist(terms_to_vars), seq_along(var_names))) == 0)

    # Check whether to compute variance between two observations
    if(type == "term_var_diff"){
      if(nrow(m) != 2)
        stop("Type 'term_var_diff' is only validt when two observations are passed")

      warning("type == 'term_var_diff' not tested")

      sds_res = array(NA_real_, dim = c(object$risk_set$d + 1, length(term_names)),
                      dimnames = list(NULL, term_names_org))

      term_dif = m[1, , drop = F] - m[2, , drop = F]
      for(i in seq_along(term_names)){
        for(j in seq_len(object$risk_set$d + 1))
          sds_res[j, i] = sqrt(diag(term_dif[, terms_to_vars[[i]], drop = F] %*%
                                      object$V_t_d_s[j, terms_to_vars[[i]], terms_to_vars[[i]]] %*%
                                      t(term_dif[, terms_to_vars[[i]], drop = F])))
      }

      return(sds_res)
    }

    # Predict terms
    terms_res = array(NA_real_, dim = c(object$risk_set$d + 1, nrow(new_data), length(term_names)),
                      dimnames = list(NULL, NULL, term_names_org))

    sds_res = if(sds) terms_res else NULL

    for(i in seq_along(term_names)){
      terms_res[, , i] = object$a_t_d_s[ ,terms_to_vars[[i]], drop = F] %*% t(m[, terms_to_vars[[i]], drop = F])
      if(!sds)
        next

      for(j in seq_len(object$risk_set$d + 1))
        sds_res[j, , i] = sqrt(diag(m[, terms_to_vars[[i]], drop = F] %*%
                                      object$V_t_d_s[j, terms_to_vars[[i]], terms_to_vars[[i]]] %*%
                                      t(m[, terms_to_vars[[i]], drop = F])))
    }

    return(list(terms = terms_res,
                sds = sds_res))
  }

  # Check if start and stop is provided. If so we need to use these
  # if not, we predict for the sample period
  if(all(c(tstart, tstop) %in% colnames(new_data))){ # TODO: Better way to test if start and stop times are provided
    # Check order of random walk
    if(object$order > 1)
      warning("Predict not test with new data for order ", object$order)

    # Find min start. Throw error if before time zero
    start = new_data[, tstart]
    stop_ = new_data[, tstop]

    if(min(start) < object$times[1])
      stop("First start time is before time zero")

    # Make prediction of covariates if last stop is beyond last period
    parems = object$a_t_d_s
    times = object$times

    max_stop = max(stop_)
    max_times = tail(object$times, 1)

    if(max_stop > tail(object$times, 1)){
      last_gab = diff(tail(object$times, 2))
      new_times = max_times + last_gab*(1:ceiling((max_stop - max_times) / last_gab))

      n_cols = dim(parems)[2]
      parems = rbind(parems, matrix(NA_real_, nrow = length(new_times), ncol = n_cols))
      if(object$order > 1)
        warning("Currently forecasting wihtout drift from higher than first order effects")
      for(t in seq_along(new_times) + length(times)){
        parems[t, ] = parems[t - 1, ]
        # parems[t, ] = object$F_ %*% parems[t - 1, ]
      }

      #       parems = c(c(t(parems)), rep(tail(parems, 1), length(new_times)))
      #       parems = matrix(parems, ncol = n_cols, byrow = T)
      times = c(times, new_times)
    }

    parems = parems[, 1:(dim(parems)[2] / object$order)] # We only need the current estimates

    # Round if needed. Throw error if so
    int_start = findInterval(start, times)
    if(any(start - times[int_start] > 0))
      warning("Some start times are rounded down")

    int_stop_ = findInterval(stop_, times + sqrt(.Machine$double.eps))
    if(any(times[int_stop_] - stop_ > 0))
      warning("Some stop times are rounded up")

    # Make function to predict for each observations
    # assume that covariates do not change
    hazard_func = object$hazard_func
    tmp_func = function(x_, sta_, sto_){
      survival_probs = 1 - sapply(sta_:sto_, function(t)
        hazard_func(parems[t, ] %*% x_))
      1 - prod(survival_probs)
    }

    # Compute hazard
    if(use_parallel){
      tryCatch({
        no_cores <- min(parallel::detectCores() - 1, nrow(m))
        cl <- parallel::makeCluster(no_cores)
        parallel::clusterExport(cl, c("parems", "hazard_func"),
                                envir = environment())

        fits = parallel::parRapply(cl = cl, data.frame(sta_ = int_start, sto_ = int_stop_, x_ = m),
                                   function(row_){
                                     tmp_func(x_ = row_[-(1:2)], sta_ = row_[1], sto_ = row_[2])
                                   })

      }, finally = { parallel::stopCluster(cl)})
    }
    else{
      fits = apply(data.frame(sta_ = int_start, sto_ = int_stop_, x_ = m), 1,
                   function(row_){
                     tmp_func(x_ = row_[-(1:2)], sta_ = row_[1], sto_ = row_[2])
                   })
    }

    return(list(fits = fits, start = times[int_start], stop = times[int_stop_]))
  }

  if(nrow(m) > 1)
    stop("Not implemented for more than one observation")

  m = m[1, ]
  var_ind = seq_along(m) # ADDED

  fits = t(object$hazard_func(object$a_t_d_s[, var_ind] %*% m)) # changed
  if(method == "backtransform"){
    sds = sqrt(apply(object$V_t_d_s, 3, function(V) diag(m %*% V[var_ind, var_ind] %*% m)))
    lbs = t(object$hazard_func(object$a_t_d_s[, var_ind] %*% m - sds * qnorm(level)))
    ubs = t(object$hazard_func(object$a_t_d_s[, var_ind] %*% m + sds * qnorm(level)))
  }
  else if(method == "delta"){
    first_deriv = apply(object$a_t_d_s[, var_ind], 1, function(b) object$hazard_first_deriv(b, m))
    sds = sapply(seq_len(nrow(object$a_t_d_s)), function(i)
      sqrt(first_deriv[, i] %*% object$V_t_d_s[var_ind, var_ind, i] %*%  first_deriv[, i]))

    lbs = fits - sds * qnorm(level)
    ubs = fits + sds * qnorm(level)
  }
  else
    stop("Method not implemented")

  list(fits = fits, lbs = lbs, ubs = ubs, times = object$times)
}
