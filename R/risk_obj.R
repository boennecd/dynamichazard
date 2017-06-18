#' Get the risk set at each bin over an equal distance grid
#' @param Y Vector of outcome variable
#' @param by Length of each bin
#' @param max_T Last observed time
#' @param id Vector with ids where entries match with outcomes \code{Y}
#' @param is_for_discrete_model \code{TRUE}/\code{FALSE} for whether the model outcome is discrete. For example, a logit model is discrete whereas what is coined an exponential model in this package is a dynamic model
#' @param n_threads Set to a value greater than one to use \code{\link{mclapply}} to find the risk object
#' @param min_chunk Minimum chunk size of ids to use when parallel version is used.
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
get_risk_obj = function(
  Y, by, max_T, id, is_for_discrete_model = T, n_threads = 1,
  min_chunk = 5000){
  if(n_threads > 1){
    unique_ids <- unique(id)
    n_ids <- length(unique_ids)
  }

  Y <- unclass(Y)
  start = Y[, 1]
  stop = Y[, 2]
  event = Y[, 3]
  is_ev <- event == 1

  start_order = order(start, method = "radix") - 1L
  max_T = if(missing(max_T)) min(max(stop[is_ev]), max(stop[-is_ev])) else max_T
  order_by_id_and_rev_start = order(id, -start, method = "radix") - 1L

  min_start <- as.double(min(start))
  storage.mode(max_T) <- "double"
  storage.mode(by) <- "double"
  event_times_in <- seq(min_start + by, max_T, by)
  tmp_n <- length(event_times_in)
  if(event_times_in[tmp_n] < max_T)
    event_times_in <- c(event_times_in, event_times_in[tmp_n] + by)

  if(n_threads == 1 || min_chunk > n_ids){
    return(
      get_risk_obj_rcpp(
        start = start, stop = stop, event = event,
        by = by, start_order = start_order,
        max_T = max_T,
        order_by_id_and_rev_start = order_by_id_and_rev_start, id = id,
        is_for_discrete_model = is_for_discrete_model,
        min_start = min_start, event_times_in = event_times_in)
    )
  }

  # Find number of tasks
  n_tasks <- min(ceiling(n_ids / min_chunk), 4 * n_threads)
  tasks <- split(unique_ids, cut(seq_along(unique_ids), n_tasks, labels = FALSE))

  # Find subset of risk sets
  out <- parallel::mclapply(tasks, function(ids) {
    my_indx <- which(id %in% ids) - 1L
    my_start_order <- intersect(start_order, my_indx)
    my_order_by_id_and_rev_start <- intersect(order_by_id_and_rev_start,  my_indx)

    local_res <- get_risk_obj_rcpp(
      start = start, stop = stop, event = event,
      by = by, start_order = my_start_order,
      max_T = max_T,
      order_by_id_and_rev_start = my_order_by_id_and_rev_start, id = id,
      is_for_discrete_model = is_for_discrete_model,
      min_start = min_start, event_times_in = event_times_in)

    local_res$is_event_in[-(my_order_by_id_and_rev_start + 1)] <- NA_integer_
    local_res
    })

  # Combine results
  final <- out[[1]]
  final$risk_sets <-
    lapply(seq_along(event_times_in), function(i){
      do.call(c, unname(lapply(out, function(x){
        x$risk_sets[[i]]
      })))
    })

  for(i in 2:length(out)){
    is_new <- which(!is.na(out[[i]]$is_event_in))
    final$is_event_in[is_new] <- out[[i]]$is_event_in[is_new]
  }

  final
}

######
# Get expression to permutate data and risk set

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
# Get expression to sort data and risk set by start and stop time

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
