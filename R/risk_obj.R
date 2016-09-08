#' Returns the discretized risk set
#' @export
get_risk_obj = function(Y, by, max_T, id, is_for_discrete_model = T)
  get_risk_obj_rcpp(
    start = Y[, 1], stop = Y[, 2], event = Y[, 3],
    by = by, start_order = order(Y[, 1]) - 1,
    max_T = ifelse(missing(max_T), max(Y[Y[, 3] == 1, 2]), max_T),
    order_by_id_and_rev_start = order(id, -Y[, 1]) - 1, id = id,
    is_for_discrete_model = is_for_discrete_model)
