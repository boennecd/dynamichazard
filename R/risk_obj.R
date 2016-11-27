#' Get the risk set at each bin over an equal distance grid
#' @param Y outcome variable
#' @param by Length of each bin
#' @param max_T Last observed time
#' @param id Vector with ids where entries match with outcomes \code{Y}
#' @param is_for_discrete_model \code{TRUE}/\code{FALSE} for whether the model outcome is discrete. For example, a logit model is discrete whereas what is coined an exponential model in this package is a dynamic model
#'
#' @return
#' A list with the following elements:
#' \tabular{ll}{
#' \code{risk_sets}\verb{ }\tab List of lists with one for each bin. Each of the sub list have indicies that coresponds to the entries of \code{Y} that are at risk in the bin \cr
#' \code{min_start}\verb{ }\tab Start time of the first bin \cr
#' \code{I_len}\verb{ }\tab Length of each bin \cr
#' \code{d}\verb{ }\tab Number of bins \cr
#' \code{is_event_in}\verb{ }\tab Indicies for which bin an observation \code{Y} is an event. \code{-1} if the individual does not die in any of the bins \cr
#' \code{is_for_discrete_model}\verb{ }\tab Value of \code{is_for_discrete_model} argument
#' }
#' @export
get_risk_obj = function(Y, by, max_T, id, is_for_discrete_model = T)
  get_risk_obj_rcpp(
    start = Y[, 1], stop = Y[, 2], event = Y[, 3],
    by = by, start_order = order(Y[, 1]) - 1,
    max_T = ifelse(missing(max_T), max(Y[Y[, 3] == 1, 2]), max_T),
    order_by_id_and_rev_start = order(id, -Y[, 1]) - 1, id = id,
    is_for_discrete_model = is_for_discrete_model)
