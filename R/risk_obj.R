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
get_risk_obj = function(Y, by, max_T, id, is_for_discrete_model = T, n_threads = 1){
  if(n_threads > 1){
    min_ids <- 2000
    unique_ids <- unique(id)
    n_ids <- length(unique_ids)
  }

  start = Y[, 1]
  stop = Y[, 2]
  event = Y[, 3]

  start_order = order(Y[, 1]) - 1L
  max_T = ifelse(missing(max_T), max(Y[Y[, 3] == 1, 2]), max_T)
  order_by_id_and_rev_start = order(id, -Y[, 1]) - 1L

  min_start <- as.double(min(start))
  storage.mode(max_T) <- "double"
  storage.mode(by) <- "double"
  event_times_in <- seq(min_start + by, max_T, by)
  if(tail(event_times_in, 1) < max_T)
    event_times_in <- c(event_times_in, tail(event_times_in, 1) + by)

  if(n_threads == 1 || min_ids > n_ids){
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
  n_tasks <- ceiling(n_ids / min_ids)
  tasks <- split(unique_ids, cut(seq_along(unique_ids), n_tasks, labels = FALSE))

  # Find subset of risk sets
  # cl <- parallel::makeCluster(min(n_threads, n_tasks), outfile = "tmp.txt") # TODO: remove outfile
  # on.exit(parallel::stopCluster(cl))
  #
  # parallel::clusterExport(cl, varlist = list(
  #   "id", "start_order", "order_by_id_and_rev_start",
  #   "start", "stop", "event", "by", "max_T", "min_start",
  #   "event_times_in", "is_for_discrete_model"),
  #   envir = environment())

  out <-  list()
  for(i in seq_along(tasks)){
    print(i)
    ids <- tasks[[i]]

    my_indx <- which(id %in% ids) - 1L
    my_start_order <- intersect(start_order, my_indx)
    my_order_by_id_and_rev_start <- intersect(order_by_id_and_rev_start,  my_indx)


    # TODO: remove
    stopifnot(all(id[my_start_order +1] %in% ids))
    stopifnot(!any(id[-(my_start_order +1)] %in% ids))

    stopifnot(all(id[my_order_by_id_and_rev_start + 1] %in% ids))
    stopifnot(!any(id[-(my_order_by_id_and_rev_start +1)] %in% ids))
    print(c(length(start), length(stop), length(id), length(ids)))
    if(length(my_start_order) != length(my_order_by_id_and_rev_start))
      print("Boh")

    out[[i]] <- get_risk_obj_rcpp(
      start = start, stop = stop, event = event,
      by = by, start_order = my_start_order,
      max_T = max_T,
      order_by_id_and_rev_start = my_order_by_id_and_rev_start, id = id,
      is_for_discrete_model = is_for_discrete_model,
      min_start = min_start, event_times_in = event_times_in)

    gc()
  }

  out

  # out <- parallel::parLapply(cl, tasks, function(ids){
  #   print("HEY") # TODO: remove
  #   print(ids) # TODO: remove
  #   my_indx <- which(id %in% ids) - 1L
  #   my_start_order <- intersect(start_order, my_indx)
  #   my_order_by_id_and_rev_start <- intersect(order_by_id_and_rev_start,  my_indx)
  #
  #   get_risk_obj_rcpp(
  #     start = start, stop = stop, event = event,
  #     by = by, start_order = my_start_order,
  #     max_T = max_T,
  #     order_by_id_and_rev_start = my_order_by_id_and_rev_start, id = id,
  #     is_for_discrete_model = is_for_discrete_model,
  #     min_start = min_start, event_times_in = event_times_in)
  # })
}

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
