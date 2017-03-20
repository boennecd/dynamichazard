#' Get the risk set at each bin over an equal distance grid
#' @param Y Vector of outcome variable
#' @param by Length of each bin
#' @param max_T Last observed time
#' @param id Vector with ids where entries match with outcomes \code{Y}
#' @param is_for_discrete_model \code{TRUE}/\code{FALSE} for whether the model outcome is discrete. For example, a logit model is discrete whereas what is coined an exponential model in this package is a dynamic model
#'
#' @return
#' A list with the following elements:
#' \describe{
#' \item{\code{risk_sets}}{List of lists with one for each bin. Each of the sub lists have indices that corresponds to the entries of \code{Y} that are at risk in the bin }
#' \item{\code{min_start}}{ Start time of the first bin }
#' \item{\code{I_len}}{ Length of each bin }
#' \item{\code{d}}{ Number of bins }
#' \item{\code{is_event_in}}{ Indices for which bin an observation \code{Y} is an event. \code{-1} if the individual does not die in any of the bins }
#' \item{\code{is_for_discrete_model}}{ Value of \code{is_for_discrete_model} argument}
#' }
#' @export
get_risk_obj = function(Y, by, max_T, id, is_for_discrete_model = T)
  get_risk_obj_rcpp(
    start = Y[, 1], stop = Y[, 2], event = Y[, 3],
    by = by, start_order = order(Y[, 1]) - 1,
    max_T = ifelse(missing(max_T), max(Y[Y[, 3] == 1, 2]), max_T),
    order_by_id_and_rev_start = order(id, -Y[, 1]) - 1, id = id,
    is_for_discrete_model = is_for_discrete_model)

get_permu_data_exp <- function(data, risk_obj, weights){
  data <- deparse(substitute(data))
  risk_obj <- deparse(substitute(risk_obj))
  weights <- deparse(substitute(weights))

  txt_exp <-
    "permu <- sample(nrow(data[[1]]), replace = F)

    for(i in seq_along(data))
      data[[i]] <- data[[i]][permu, , drop = F]

    for(i in seq_along(risk_obj$risk_sets))
      risk_obj$risk_sets[[i]] <-
        structure(sort(match(risk_obj$risk_sets[[i]], permu)), org = risk_obj$risk_sets[[i]])

    risk_obj$is_event_in <- risk_obj$is_event_in[permu]

    weights <- weights[permu]

    rm(i)"

  txt_exp <- gsub("data", data, txt_exp)
  txt_exp <- gsub("risk_obj", risk_obj, txt_exp)
  txt_exp <- gsub("weights", weights, txt_exp)

  parse(text = txt_exp)
}

get_permu_data_rev_exp <- function(data, risk_obj, weights){
  data <- deparse(substitute(data))
  risk_obj <- deparse(substitute(risk_obj))
  weights <- deparse(substitute(weights))

  txt_exp <-
    "org_order <- order(permu)

    for(i in seq_along(data))
      data[[i]] <- data[[i]][org_order, , drop = F]

    for(i in seq_along(risk_obj$risk_sets))
      risk_obj$risk_sets[[i]] <- attr(risk_obj$risk_sets[[i]], 'org')

    risk_obj$is_event_in <- risk_obj$is_event_in[org_order]

    weights <- weights[org_order]

    rm(permu, org_order, i)"

  txt_exp <- gsub("data", data, txt_exp)
  txt_exp <- gsub("risk_obj", risk_obj, txt_exp)
  txt_exp <- gsub("weights", weights, txt_exp)

  parse(text = txt_exp)
}
