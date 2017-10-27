#' @title Summarizing Dynamic Hazard Models Fits
#'
#' @param object object returned from \code{\link{ddhazard}}.
#' @param x object returned from \code{summary.ddhazard}.
#' @param var_indices variable indices to print for time-varying effects.
#' @param max_print maximum number of time points to print coefficients at.
#' @param digits number of digits to print.
#' @param ... not used.
#'
#' @description
#' The \code{sd} printed for time-varying effects are point-wise standard deviations from the smoothed covariance matrices.
#'
#' @export
summary.ddhazard <- function(
  object, var_indices = 1:ncol(object$state_vecs), max_print = 10, ...){
  ans <- object[c("call", "model", "method", "n_iter")]

  state_vecs <- object$state_vecs
  state_vars <- object$state_vars
  state_vars <-
    if(ncol(object$state_vecs) > 1)
      t(apply(state_vars, 3, diag)) else if(ncol(object$state_vecs) == 1)
        as.matrix(apply(state_vars, 3, diag), ncol = 1) else
          state_vecs

  if(length(state_vecs) > 0 && length(var_indices) > 0 && max_print > 0){
    time_indices <- if(nrow(state_vecs) > max_print){
      max_print <- max_print - 1
      i1 <- (1:nrow(state_vecs) - 1) * max_print / nrow(state_vecs)
      i2 <- tapply(i1, as.integer(i1), max)
      c(1, which(i1 %in% i2))
    } else
      1:nrow(state_vecs)
    out <- cbind(state_vecs[time_indices, var_indices, drop = F],
                 sqrt(state_vars[time_indices, var_indices, drop = F]))

    colnames(out) <- c(colnames(out)[seq_along(var_indices)],
                       rep("  sd ", length(var_indices)))

    out <- out[, c(sapply(seq_along(var_indices), rep, times = 2)) +
                 rep(c(0, length(var_indices)), length(var_indices))]

    rownames(out) <- format(object$times[time_indices], justify = "left")

    ans$coefficients <- out
  }

  ans$Q <- object$Q
  ans$n_id <- length(unique(object$id))
  ans$n_events <- if(!is.null(object$risk_set))
    sum(object$risk_set$is_event_in > -1) else NULL
  ans$fixed_effects <- object$fixed_effects
  class(ans) <- "summary.ddhazard"
  ans
}
