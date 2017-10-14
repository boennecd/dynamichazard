#' @title Hat values for \code{\link{ddhazard}}
#'
#' @description Computes hat-"like" values from usual L2 penalized binary regression
#'
#' @param model A fit from \code{\link{ddhazard}}
#' @param ... Not used
#'
#' @details
#' Computes hat-"like" values in each interval for each individual at risk in the interval. See the ddhazard vignette for details
#'
#' @return
#' A list of matrices. Each matrix has three columns: the hat values, the row number of the original data point and the id the row belongs to
#'
#' @seealso
#' \code{\link{ddhazard}}
#'
#' @export
hatvalues.fahrmeier_94 <- function(model, ...){
  if(!model$model %in% c("logit"))
    stop("Functions for model '",  model$model, "' is not implemented")

  if(is.null(model) || is.null(model$risk_set))
    stop("Missing risk set or data to compute residuals. Call ddhazard again with Control values set to save these")

  # Get design matrix, save the start and stop times and include the fixed
  # effects
  design_mat <- get_design_matrix(
    data = model$data, Terms = model$terms, xlev = model$xlev,
    has_fixed_intercept = model$has_fixed_intercept)
  tstart <- design_mat$Y[, 1]
  tstop <- design_mat$Y[, 2]
  design_mat = cbind(design_mat$X, design_mat$fixed_terms)

  # Get the coeffecient matrix including fixed effects
  coefs <- model$state_vecs
  n_varying <- length(coefs)
  if(n_varying > 0){
    coefs = coefs[-1, ] # remove first row it is the initial state space vector
    n_varying <- dim(coefs)[2] / model$order
    coefs = coefs[, 1:n_varying] # We only need the current estimates (relevant for higher than 1. order)
  }
  if(length(model$fixed_effects) > 0)
    coefs <- cbind(coefs, sapply(model$fixed_effects, rep, times = nrow(coefs)))

  # Compute hat values
  V_inv <- matrix(0.0, nrow = ncol(coefs), ncol = ncol(coefs))
  res <- list()
  for(i in seq_along(model$risk_set$risk_sets)){
    start_ = model$times[i]
    stop_ = model$times[i + 1]
    r_set = model$risk_set[[1]][[i]]

    # Find design matrix and coeffecients
    design_mat_sub <- design_mat[r_set, ]
    t_design_mat_sub <- t(design_mat_sub)
    coefs_cur <- coefs[i, ]

    if(n_varying > 0)
      V_inv[1:n_varying, 1:n_varying] <-
        solve(model$state_vars[1:n_varying, 1:n_varying, i + 1])

    # Find start and stop times + the linear predictors
    tsta = pmax(tstart[r_set], start_)
    tsto <- pmin(tstop[r_set], stop_)
    etas = design_mat_sub %*% coefs_cur

    # Compute the variances
    vars <- mapply(model$var_func, eta = etas, tstart = tsta, tstop = tsto) + model$control$denom_term

    # Compute the hat values
    # H <- solve(t_design_mat_sub %*% diag(vars) %*% design_mat_sub  + V_inv)
    H <- solve(t_design_mat_sub %*% (vars * design_mat_sub) + V_inv)

    # library(microbenchmark)
    # microbenchmark(
    #   t_design_mat_sub %*% diag(vars) %*% design_mat_sub,
    #   t_design_mat_sub %*% (vars * design_mat_sub))

    # H <- diag(design_mat_sub %*% H %*% t_design_mat_sub) * vars
    H <- colSums(t_design_mat_sub * (H %*% t_design_mat_sub)) * vars

    # library(microbenchmark)
    # microbenchmark(
    #   diag(design_mat_sub %*% H %*% t_design_mat_sub) * vars,
    #   colSums(t_design_mat_sub * (H %*% t_design_mat_sub)) * vars)

    tmp <- cbind("hat_value" = H,
                 row_num = r_set,
                 id = model$id[r_set])
    row.names(tmp) <- NULL
    res[[i]] <- tmp
  }

  res
}

