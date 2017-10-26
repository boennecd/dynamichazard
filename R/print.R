# from boot.print
#' @title Summary statistics for a \code{ddhazard_boot} object
#'
#' @param x returned object from a \code{\link{ddhazard_boot}} call.
#' @param digits the number of digits to be printed in the summary statistics.
#' @param index indices indicating for which elements of the bootstrap output summary statistics are required.
#' @param ... not used.
#'
#' @description
#' Arguments have the same effects as for an object from a \code{\link{boot}} call. See \code{\link[=print.boot]{print}}.
#'
#' @seealso
#' \code{\link{ddhazard_boot}}
#'
#' @export
print.ddhazard_boot <-
  function (x, digits = getOption("digits"), index = 1L:ncol(boot.out$t), ...)
  {
    boot.out <- x
    sim <- boot.out$sim
    cl <- boot.out$call
    t <- matrix(boot.out$t[, index], nrow = nrow(boot.out$t))
    allNA <- apply(t, 2L, function(t) all(is.na(t)))
    ind1 <- index[allNA]
    index <- index[!allNA]
    t <- matrix(t[, !allNA], nrow = nrow(t))
    rn <- boot.out$t_names
    if (length(index) == 0L)
      op <- NULL
    else {
      t0 <- boot.out$t0[index, drop = F]
      op <- cbind(
        t0,
        apply(t, 2L, mean, na.rm = TRUE) - t0,
        apply(t, 2L, mean, na.rm = TRUE, trim = .025) - t0,
        sqrt(apply(t, 2L, function(t.st) var(t.st[!is.na(t.st)]))))
      dimnames(op) <-
        list(rn, c("original", " bias  "," bias (truncated)"," std. error"))
    }

    cat("Bootstrap Statistics :\n")
    if (!is.null(op))
      print(op, digits = digits)
    invisible(boot.out)
  }


#' @title Print function for \code{ddhazard} result
#'
#' @param x object returned from \code{\link{ddhazard}}.
#' @param var_indices variable indices to print for time-varying effects.
#' @param time_indices time intervals to print for time-varying effects.
#' @param digits number of digits to print.
#' @param ... not used.
#'
#' @description
#' The \code{sd} printed for time-varying effects are point-wise standard deviations from the smoothed covariance matrices.
#'
#' @export
print.ddhazard<- function(
  x, var_indices = 1:ncol(x$state_vecs), time_indices = 1:nrow(x$state_vecs),
  digits = getOption("digits"), ...){
  cat("Formula:\n", deparse(x$formula), "\n", sep = "")

  cat("\nEstimated with ", x$method, " in ", x$n_iter, " iterations of the EM algorithm\n",
      sep = "")

  state_vecs <- x$state_vecs
  state_vars <- x$state_vars
  state_vars <-
    if(ncol(x$state_vecs) > 1)
      t(apply(state_vars, 3, diag)) else if(ncol(x$state_vecs) == 1)
        as.matrix(apply(state_vars, 3, diag), ncol = 1) else
          state_vecs

  if(length(state_vecs) > 0 && length(var_indices) > 0 &&
      length(time_indices) > 0){
    out <- cbind(state_vecs[time_indices, var_indices, drop = F],
                 sqrt(state_vars[time_indices, var_indices, drop = F]))

    colnames(out) <- c(colnames(out)[seq_along(var_indices)],
                       rep("  sd ", length(var_indices)))

    out <- out[, c(sapply(seq_along(var_indices), rep, times = 2)) +
                 rep(c(0, length(var_indices)), length(var_indices))]

    rownames(out) <- paste0("t", 1:nrow(out) - 1)

    cat("\nEstimated time-varying effects and point-wise standard deviation:\n")
    print(out, digits = digits)
  }

  fixed_effects <- x$fixed_effects
  if(length(fixed_effects) > 0){
    cat("\nFixed effects are estimated in the ", x$control$fixed_terms_method,
        ". The estimates are:\n", sep = "")
    print(fixed_effects, digits = digits)
  }

  invisible(x)
}
