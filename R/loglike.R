#' @title Log likelihood of smoothed state vector of a \code{fahrmeier_94} object
#' @description
#' Computes the log likelihood of (a potentially new) data set given the estimated:
#' \deqn{E_{\theta}(\alpha_1 | y_{1:d}), E_{\theta}(\alpha_{2} | y_{1:d}), ..., E_{\theta}(\alpha_{d} | y_{1:d})}
#'
#' from the \code{fahrmeier_94} class object. Note that this is not the log likelihood of the observed data given the outcome.
#'
#' @param object an object of class \code{fahrmeier_94}.
#' @param data new data to evaluate the likelihood for.
#' @param id the individual identifiers as in \code{\link{ddhazard}}.
#' @param ... unused.
#'
#' @export
logLik.fahrmeier_94 = function(object, data = NULL, id, ...){
  data <- if(!is.null(object$data)) object$data else data
  if(is.null(data))
    stop("data is needed to compute log likelihood. Please, pass the data set used in 'ddhazard' call")

  X <- get_design_matrix(
    data = data, Terms = object$terms, xlev = object$xlev,
    has_fixed_intercept = object$has_fixed_intercept)
  X$X <- t(X$X)

  fixed_effects_offsets <- if(ncol(X$fixed_terms) == 0)
    rep(0, nrow(X$fixed_terms)) else
      X$fixed_terms %*% object$fixed_effects

  risk_obj <- object$risk_set
  if(is.null(risk_obj)){
    if(missing(id))
      stop("id need to compute log likelihood. Please, pass the id used in 'ddhazard' call")

    if(object$model %in% exp_model_names){
      is_for_discrete_model <- F

    } else if (object$model == "logit"){
      is_for_discrete_model <- T

    } else
      stop("logLik not implemented for model '", object$model, "'")

    risk_obj <- get_risk_obj(Y = X$Y, by = unique(diff(object$times)),
                             max_T = max(object$times), is_for_discrete_model = is_for_discrete_model,
                             id = id)
  }

  val <- logLike_cpp(X = X$X, risk_obj = risk_obj, F = object$F_,
                     Q_0 = object$Q_0, Q = object$Q, a_t_d_s = t(object$state_vecs),
                     tstart = X$Y[, 1], tstop = X$Y[, 2], order_ = object$order,
                     model = object$model, fixed_effects_offsets = fixed_effects_offsets)

  attr(val, "prior_loglike") <- val[2]
  val <- val[1]

  if(object$est_Q_0)
    warning("parameters for Q_0 are not included in attribute df")

  # n_parems <- ncol(object$state_vecs) / object$order
  # attr(val, "df") <- n_parems * object$order +  # from a_0
  #   n_parems * (n_parems + 1) / 2 +  # from Q
  #   length(object$fixed_effects) # from fixed effects
  class(val) <- "logLik"

  val
}

#' @title Log-Likelihood of a \code{PF_clouds} object
#' @param object an object of class \code{PF_clouds}
#' @param ... unused.
#' @description
#' Computes the log-likelihood using the forward filter clouds. See the particle_filter vignette for details.
#'
#' @return
#' The log-likelihood value given the observed data and set of parameter used when simulating the clouds.
#'
#' @export
logLik.PF_clouds <- function(object, ...){
  sum(tail(
    sapply(lapply(object$forward_clouds,
                  "[[", "log_unnormalized_weights"),
           function(x){
             .max <- max(x)
             log(sum(exp(x - .max))) + .max - log(length(x))
           }),
    -1))
}
