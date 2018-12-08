# from boot.print
#' @title Summary Statistics for a ddhazard_boot Object
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


#' @export
print.ddhazard<- function(x, ...){
  cat("Call:", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "\n")

  cat("\n", sQuote(x$model), " model fitted with the ", sQuote(x$method),
      " method in ", x$n_iter, " iterations of the EM algorithm.\n", sep = "")

  invisible(x)
}

#' @rdname summary.ddhazard
#' @export
print.summary.ddhazard <- function(x, digits = getOption("digits"), ...){
  old <- getOption("digits")
  on.exit(options(digits = old))
  options(digits = digits)

  print.ddhazard(x)
  cat("\n")

  if(!is.null(x$coefficients)){
    cat("Smoothed time-varying coefficients are:\n")
    print(x$coefficients)
    cat("\n")
  }

  if(!is.null(x$Q)){
    cat("The estimated diagonal entries of the covariance matrix in the state equation are:\n")
    print(diag(x$Q))
    cat("\n")
  }

  if(!is.null(x$fixed_effects) && length(x$fixed_effects) > 0){
    cat("The estimated fixed effects are:\n")
    print(x$fixed_effects)
    cat("\n")
  }

  cat(
    x$n_id, " individuals used in estimation",
    if(!is.null(x$n_events))
      paste(" with", x$n_events, "observed events") else "",
    ".\n", sep = "")

  invisible(x)
}

#' @importFrom utils tail
#' @export
print.PF_EM <- function(x, ...){
  cat("Call:\n", paste0(deparse(x$call), collapse = "\n"),
      "\n\n", sep = "")

  cat("Model estimated in ", x$n_iter, " iterations of the EM algorithm. ",
      "The log-likelihood in the last iteration is ", tail(x$log_likes, 1),
      ".\n", sep = "")

  invisible(x)
}

#' @export
print.PF_clouds <- function(x, ...){
  cat("Particle clouds with ", length(x$forward_clouds),
      " forward filter clouds ", length(x$backward_clouds),
      " backward filter clouds and ", length(x$smoothed_clouds),
      " clouds from smoothing.\n", sep = "")

  invisible(x)
}

#' @importFrom utils head
#' @export
print.ddhazard_space_errors <- function(x, ...){
  type <- if(x$standardize)
    "Standardized state space errors" else "State space errors"

  n <- 5
  cat(type, " from a model with ", nrow(x$residuals), " periods and ",
      ncol(x$residuals), " coefficients. The errors in the first ", n,
      " time periods are:\n", sep = "")
  errs <- head(x$residuals, n)

  dimnames(errs) <- list("Period" = 1:n, "Coefficient" = colnames(errs))
  print(errs)

  invisible(x)
}

#' @export
print.ddsurvcurve <- function(x, ...){
  cat("Predicted survival curve from call:\t",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "\n")
}
