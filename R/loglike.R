#' @export
logLik.fahrmeier_94 = function(object, data = NULL, id, ...){
  data <- if(!is.null(object$data)) object$data else data
  if(is.null(data))
    stop("data need to compute log likelihood. Please, pass the data set used in 'ddhazard' call")

  X <- get_design_matrix(object$formula, data)
  X$X <- t(X$X)

  risk_obj <- object$risk_set
  if(is.null(risk_obj)){
    if(missing(id))
      stop("id need to compute log likelihood. Please, pass the id used in 'ddhazard' call")

    if(object$model == "exponential"){
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
                     model = object$model)

  attr(val, "prior_loglike") <- val[2]
  val <- val[1]

  if(object$est_Q_0)
    warning("parameters for Q_0 are not included in attribute df")

  n_parems <- ncol(object$state_vecs) / object$order
  attr(val, "df") <- n_parems * object$order +  # from a_0
    n_parems * (n_parems + 1) / 2 # from Q
  class(val) <- "logLik"

  val
}
