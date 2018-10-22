#' @title Log Likelihood of Mean Path of ddhazard Object
#'
#' @description
#' Computes the log likelihood of (a potentially new) data set given the estimated:
#' \deqn{E_{\theta}(\alpha_1 | y_{1:d}), E_{\theta}(\alpha_{2} | y_{1:d}), ..., E_{\theta}(\alpha_{d} | y_{1:d})}
#'
#' of the \code{ddhazard} object. Note that this is not the log likelihood of the observed data given the outcome.
#'
#' @param object an object of class \code{ddhazard}.
#' @param data new data to evaluate the likelihood for.
#' @param id the individual identifiers as in \code{\link{ddhazard}}.
#' @param ... unused.
#'
#' @examples
#'library(dynamichazard)
#'fit <- ddhazard(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
#'  Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
#'  control = ddhazard_control(method = "GMA"))
#'logLik(fit)
#'
#' @export
logLik.ddhazard = function(object, data = NULL, id, ...){
  data <- if(!is.null(object$data)) object$data else data
  if(is.null(data))
    stop("data is needed to compute log likelihood. Please, pass the data set used in 'ddhazard' call")

  X <- get_design_matrix(
    formula = object$formula,
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

    risk_obj <- get_risk_obj(
      Y = X$Y, by = unique(diff(object$times)), id = id,
      max_T = max(object$times), is_for_discrete_model = is_for_discrete_model)
  }

  val <- logLike_cpp(
    X = X$X, risk_obj = risk_obj, F = object$F_,
    Q_0 = object$Q_0, Q = object$Q, a_t_d_s = t(object$state_vecs),
    tstart = X$Y[, 1], tstop = X$Y[, 2], order_ = object$order,
    model = object$model, fixed_effects_offsets = fixed_effects_offsets)

  attr(val, "prior_loglike") <- val[2]
  val <- val[1]

  if(object$est_Q_0)
    warning("parameters for Q_0 are not included in attribute df")

  class(val) <- "logLik"
  val
}

#' @title Approximate Log-Likelihood from a Particle Filter
#' @param object an object of class \code{PF_clouds} or \code{PF_EM}.
#' @param df degrees of freedom used in the model.
#' @param nobs integer with number of individuals used to estimate the
#' model.
#' @param ... unused.
#' @description
#' Computes the approximate log-likelihood using the forward filter clouds. See
#' the \code{vignette("Particle_filtering", "dynamichazard")} for details.
#'
#' @return
#' The approximate log-likelihood value given the observed data and set of
#' parameter used when simulating the clouds. An attribute
#' \code{"P(y_t|y_{1:(t-1)})"} has the \eqn{P(y_t|y_{1:(t-1)})} terms.
#'
#' @export
logLik.PF_EM <- function(object, ...){
  q <- ncol(object$Q)
  type <- object$type
  df <- if(type == "VAR"){
    if(!is.null(object$psi))
      with(object, length(psi) + length(phi) + length(theta)) else
        q * (q + 1) / 2 + length(object$F)
  } else if(type == "RW")
    q * (q + 1) / 2 else
      stop("Not implemented for type ", sQuote(type))

  df <- df + length(object$fixed_effects)

  logLik(object$clouds, df = df)
}

##' @rdname logLik.PF_EM
##' @method logLik PF_clouds
##' @export
logLik.PF_clouds <- function(object, df = NA_real_, nobs = NA_integer_, ...){
  term <- sapply(lapply(object$forward_clouds,
                        "[[", "log_likelihood_term"),
                 function(x){
                   .max <- max(x)
                   log(sum(exp(x - .max))) + .max
                 })[-1]

  structure(sum(term), "P(y_t|y_{1:(t-1)})" = term,
            df = df, nobs = nobs, class = "logLik")
}
