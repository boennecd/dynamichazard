#' Function used to get design matrix and weigths for a static fit where
#' there are not duplicate rows but weights instead
#' @export
get_survival_case_Weigths_and_data = function(
  formula, data, by, max_T, id, init_weights, risk_obj,
  is_for_discrete_model = T){
  X_Y = get_design_matrix(formula, data)

  if(missing(init_weights))
    init_weights = rep(1, nrow(data))

  if(missing(id)){
    warning("You did not parse and ID argument. I do not hink this is what you want ...")
    id = 1:nrow(data)
  }

  if(any(colnames(data) == "weights")){
    warning("Column called weights will be replaced")
    data = data[, colnames(data) != "weights"]
  }

  compute_risk_obj <- missing(risk_obj) || is.null(risk_obj)

  if(!compute_risk_obj && nrow(data) <
     max(unlist(lapply(risk_obj$risk_sets, max))))
    stop("risk_obj has indicies out site of data. Likely the risk_set comes from a different data set")

  if(compute_risk_obj)
    risk_obj <- get_risk_obj(
      Y = X_Y$Y, by = by,
      max_T = ifelse(missing(max_T), max(X_Y$Y[X_Y$Y[, 3] == 1, 2]), max_T),
      id = id, is_for_discrete_model = is_for_discrete_model)
      # start = X_Y$Y[, 1], stop = X_Y$Y[, 2], event = X_Y$Y[, 3],
      # by = by, start_order = order(X_Y$Y[, 1]) - 1,
      # max_T = ifelse(missing(max_T), max(X_Y$Y[X_Y$Y[, 3] == 1, 2]), max_T),
      # order_by_id_and_rev_start = order(id, -X_Y$Y[, 1]) - 1, id = id)

  new_weights = rep(0, nrow(data))
  new_case_rows = data.frame()

  for(i in seq_along(risk_obj$risk_sets)){
    time_ = risk_obj$event_times[i]
    r_set = risk_obj$risk_sets[[i]]

    is_case = risk_obj$is_event_in[r_set] == i - 1

    new_case_rows = rbind(new_case_rows, cbind(
      Y = rep(1, sum(is_case)), data[r_set[is_case], ], weights = init_weights[r_set[is_case]]))

    new_weights[r_set[!is_case]] = new_weights[r_set[!is_case]] + init_weights[r_set[!is_case]]
  }

  X = cbind(Y = rep(0, nrow(data)), data, weights = new_weights)[new_weights > 0, ]
  X = rbind(X, new_case_rows)

  X
}


#' Function to make a static glm fit from a input with a \code{Surv} object
#' as the right hand site of forumla
#' @export
static_glm = function(form, data, by, max_T, id, family = "binomial", model = F, weights, risk_obj = NULL, ...){
  if(family == "binomial"){
    family <- binomial()
  } else
    stop("family '", family, "' not implemented in static_glm")

  X = get_survival_case_Weigths_and_data(
    formula = form, data = data, by = by, max_T = max_T, id = id,
    init_weights = weights, risk_obj = risk_obj)

  weights = X$weights

  form <- formula(terms(form, data = data))
  form = update(form, Y ~ ., data = X)

  glm(form, X, family = family, model = model, weights = weights, ...)
}
