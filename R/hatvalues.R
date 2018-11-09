#' @title Hat Values for ddhazard Object
#'
#' @description Computes hat-"like" values from usual L2 penalized binary regression.
#'
#' @param model a fit from \code{\link{ddhazard}}.
#' @param ... not used.
#'
#' @details
#' Computes hat-"like" values in each interval for each individual at risk in the interval. See the \code{vignette("ddhazard", "dynamichazard")} vignette for details.
#'
#' @return
#' A list of matrices. Each matrix has three columns: the hat values, the row number of the original data point and the id the row belongs to.
#'
#' @seealso
#' \code{\link{ddhazard}}
#'
#' @examples
#'library(dynamichazard)
#'fit <- ddhazard(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3000,
#'  Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 100,
#'  control = ddhazard_control(method = "GMA"))
#'hvs <- hatvalues(fit)
#'head(hvs[[1]])
#'head(hvs[[2]])
#'
#' @export
hatvalues.ddhazard <- function(model, ...){
  if(!model$model %in% c("logit"))
    stop("Functions for model ",  sQuote(model$model), " is not implemented")

  if(is.null(model) || is.null(model$risk_set))
    stop("Missing risk set or data to compute residuals. Call ddhazard again with Control values set to save these")

  # get design matrix, save the start and stop times and include the fixed
  # effects
  design_mat <- get_design_matrix(
    formula = model$formula,
    data = model$data, Terms = model$terms, xlev = model$xlev,
    has_fixed_intercept = model$has_fixed_intercept)
  tstart <- design_mat$Y[, 1]
  tstop <- design_mat$Y[, 2]
  design_mat = cbind(design_mat$X, design_mat$fixed_terms)

  # get the coeffecient matrix including fixed effects
  coefs <- model$state_vecs
  n_varying <- length(coefs)
  if(n_varying > 0){
    coefs = coefs[-1, ] # remove the initial state space vector
    n_varying <- dim(coefs)[2] / model$order
    coefs = coefs[, 1:n_varying] # we only need the current estimates --
                                 # relevant for higher than 1. order RW

  }

  if(length(model$fixed_effects) > 0)
    coefs <- cbind(coefs, sapply(model$fixed_effects, rep, times = nrow(coefs)))

  # compute hat values
  V_inv <- matrix(0.0, nrow = ncol(coefs), ncol = ncol(coefs))
  res <- list()
  for(i in seq_along(model$risk_set$risk_sets)){
    start_ = model$times[i]
    stop_ = model$times[i + 1]
    r_set = model$risk_set[[1]][[i]]

    # find design matrix and coeffecients
    design_mat_sub <- design_mat[r_set, ]
    t_design_mat_sub <- t(design_mat_sub)
    coefs_cur <- coefs[i, ]

    if(n_varying > 0)
      V_inv[1:n_varying, 1:n_varying] <-
        solve(model$state_vars[1:n_varying, 1:n_varying, i + 1])

    # find start and stop times & the linear predictors
    tsta = pmax(tstart[r_set], start_)
    tsto <- pmin(tstop[r_set], stop_)
    etas = design_mat_sub %*% coefs_cur

    # compute the working weights
    V <- local({
      var <-
        model$family$var(eta = etas, at_risk_length = tsto - tsta) +
        model$control$denom_term
      mu_eta <-
        model$family$mu_eta(eta = etas, at_risk_length = tsto - tsta) +
        model$control$denom_term

      model$weights[r_set] * mu_eta^2 / var
    })

    # compute the hat values
    H <- t_design_mat_sub %*% (V * design_mat_sub) + V_inv
    H <- colSums(t_design_mat_sub * solve(H, t_design_mat_sub)) * V

    tmp <- cbind("hat_value" = H, row_num = r_set, id = model$id[r_set])
    row.names(tmp) <- NULL
    res[[i]] <- tmp
  }

  res
}
