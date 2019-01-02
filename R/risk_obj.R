#' @title Risk Set on an Equidistant Distant Grid
#' @description Get the risk set at each bin over an equidistant distant grid.
#'
#' @param Y vector of outcome variable returned from \code{\link{Surv}}.
#' @param by length of each bin.
#' @param max_T last observed time.
#' @param id vector with ids where entries match with outcomes \code{Y}.
#' @param is_for_discrete_model \code{TRUE} if the model outcome is discrete. For example, a logit model is discrete whereas what is is referred to as the exponential model in this package is a dynamic model.
#' @param n_threads set to a value greater than one to use \code{\link{mclapply}} to find the risk object.
#' @param min_chunk minimum chunk size of ids to use when parallel version is used.
#'
#' @return
#' a list with the following elements
#' \item{risk_sets}{list of lists with one for each bin. Each of the sub lists have indices that corresponds to the entries of \code{Y} that are at risk in the bin.}
#' \item{min_start}{start time of the first bin.}
#' \item{I_len}{length of each bin.}
#' \item{d}{number of bins.}
#' \item{is_event_in}{indices for which bin an observation \code{Y} is an event. \code{-1} if the individual does not die in any of the bins.}
#' \item{is_for_discrete_model}{value of \code{is_for_discrete_model} argument.}
#'
#' @examples
#'# small toy example with time-varying covariates
#'dat <- data.frame(
#'  id     = c(1, 1, 2, 2),
#'  tstart = c(0, 4, 0, 2),
#'  tstop  = c(4, 6, 2, 4),
#'  event  = c(0, 1, 0, 0))
#'
#'with(dat, get_risk_obj(Surv(tstart, tstop, event), by = 1, max_T = 6, id = id))
#'
#' @export
get_risk_obj = function(
  Y, by, max_T, id, is_for_discrete_model = T, n_threads = 1,
  min_chunk = 5000){
  Y <- unclass(Y)
  start = Y[, 1]
  stop = Y[, 2]
  event = Y[, 3]
  is_ev <- event == 1

  max_T = if(missing(max_T))
    min(max(stop[is_ev]), max(stop[!is_ev])) else max_T
  order_by_id_and_rev_start = order(id, -start, method = "radix") - 1L

  min_start <- as.double(min(start))
  storage.mode(max_T) <- "double"
  storage.mode(by) <- "double"
  event_times_in <- seq(min_start + by, max_T, by)
  tmp_n <- length(event_times_in)
  if(!isTRUE(all.equal(event_times_in[tmp_n], max_T)))
    event_times_in <- c(event_times_in, event_times_in[tmp_n] + by)

  # set exactly to boundaries where values are very close and difference is
  # likely due to floating point operations
  start_order = order(start, method = "radix") - 1L
  start <- round_if_almost_eq(start, start_order, c(min_start, event_times_in))
  stop_order = order(stop, method = "radix") - 1L
  stop <- round_if_almost_eq(stop, stop_order, c(min_start, event_times_in))

  if(n_threads > 1)
    warning(sQuote("n_threads"), " greater than one is no longer supported")

  get_risk_obj_rcpp(
    start = start, stop = stop, event = event,
    by = by, start_order = start_order,
    max_T = max_T,
    order_by_id_and_rev_start = order_by_id_and_rev_start, id = id,
    is_for_discrete_model = is_for_discrete_model,
    min_start = min_start, event_times_in = event_times_in)
}

######
# get expression to permutate data and risk set

permu_txt <-
  "for(i in seq_along(data))
      data[[i]] <- data[[i]][permu, , drop = F]

    for(i in seq_along(risk_obj$risk_sets))
      risk_obj$risk_sets[[i]] <-
        structure(sort(match(risk_obj$risk_sets[[i]], permu)), org = risk_obj$risk_sets[[i]])

    risk_obj$is_event_in <- risk_obj$is_event_in[permu]

    weights <- weights[permu]

    rm(i)"

get_permu_data_exp <- function(data, risk_obj, weights){
  data <- deparse(substitute(data))
  risk_obj <- deparse(substitute(risk_obj))
  weights <- deparse(substitute(weights))

  txt_exp <- paste0(
    "permu <- sample(nrow(data[[1]]), replace = F)

    ", permu_txt)

  txt_exp <- gsub("data", data, txt_exp)
  txt_exp <- gsub("risk_obj", risk_obj, txt_exp)
  txt_exp <- gsub("weights", weights, txt_exp)

  parse(text = txt_exp)
}

permu_rev_txt <-
  "org_order <- order(permu)

  for(i in seq_along(data))
    data[[i]] <- data[[i]][org_order, , drop = F]

  for(i in seq_along(risk_obj$risk_sets))
    risk_obj$risk_sets[[i]] <- attr(risk_obj$risk_sets[[i]], 'org')

  risk_obj$is_event_in <- risk_obj$is_event_in[org_order]

  weights <- weights[org_order]

  rm(permu, org_order, i)"

get_permu_data_rev_exp <- function(data, risk_obj, weights){
  data <- deparse(substitute(data))
  risk_obj <- deparse(substitute(risk_obj))
  weights <- deparse(substitute(weights))

  txt_exp <- permu_rev_txt

  txt_exp <- gsub("data", data, txt_exp)
  txt_exp <- gsub("risk_obj", risk_obj, txt_exp)
  txt_exp <- gsub("weights", weights, txt_exp)

  parse(text = txt_exp)
}

######
# get expression to sort data and risk set by start and stop time

get_order_data_exp <- function(data, risk_obj, weights){
  data <- deparse(substitute(data))
  risk_obj <- deparse(substitute(risk_obj))
  weights <- deparse(substitute(weights))

  txt_exp <- paste0(
    "permu <- order(X_Y[['Y']][, 2], X_Y[['Y']][, 1])

    ", permu_txt)

  txt_exp <- gsub("data", data, txt_exp)
  txt_exp <- gsub("risk_obj", risk_obj, txt_exp)
  txt_exp <- gsub("weights", weights, txt_exp)

  parse(text = txt_exp)
}

get_order_data_rev_exp <- get_permu_data_rev_exp
