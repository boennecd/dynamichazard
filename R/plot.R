#' @title Plots for \code{\link{ddhazard}}
#' @description Plot to illustrate the estimate state space variables from a \code{\link{ddhazard}} fit
#'
#' @param x Result of \code{\link{ddhazard}} call
#' @param type Type of plot. Currently, only \code{"cov"} is available for plot of the state space parameters
#' @param plot_type The \code{type} argument passed to \code{plot}
#' @param cov_index The index (indices) of the state space parameter(s) to plot
#' @param add \code{FALSE} if you want to make a new plot
#' @param xlab,ylab,ylim,col Arguments to overide defaults set in the function
#' @param ... Arguments passed to \code{plot} or \code{lines} depending on the value of \code{add}
#'
#' @details
#' Creates a plot of state variables or adds state variables to a plot with indices \code{cov_index}. Pointwise 1.96 std. confidence intervals are provided with the smoothed co-variance matrices from the fit
#'
#' @export
plot.fahrmeier_94 = function(x, xlab = "Time",
                             ylab = "Hazard",
                             type = "cov", plot_type = "l", cov_index, ylim,
                             col = "black", add = F, ...){
  if(!x$model %in% c("logit", "exponential"))
    stop("Functions for model '", x$model, "' is not implemented")

  if(type == "cov"){
    if(missing(cov_index))
      stop("Need cov index when type is equal to ", type)

    if(missing(ylab))
      ylab = colnames(x$state_vecs)[cov_index]

    lb = x$state_vecs[, cov_index] - 1.96 * sqrt(x$state_vars[cov_index, cov_index, ])
    ub = x$state_vecs[, cov_index] + 1.96 * sqrt(x$state_vars[cov_index, cov_index, ])

    if(missing(ylim))
      ylim = range(lb, ub)

    if(!add){
      plot(x$times, x$state_vecs[, cov_index], type = plot_type,
           ylim = ylim, xlab = xlab, ylab = ylab, col = col, ...)

    } else {
      lines(x$times, x$state_vecs[, cov_index], col = col, ...)

    }

    lines(x$times, lb, lty = 2, col = col)
    lines(x$times, ub, lty = 2, col = col)

    return(invisible())
  }

  stop("Type '", type, "' is not implemented for plot.fahrmeier_94")
}

#' @title State space error plot
#' @description Plot function for state space errors from \code{\link{ddhazard}} fit
#'
#' @param x Result of \code{\link[=residuals.fahrmeier_94]{residuals}} for state space errors
#' @param mod The \code{\link{ddhazard}} result used in the \code{\link[=residuals.fahrmeier_94]{residuals}} call
#' @param p_cex \code{cex} argument for the points
#' @param cov_index The indices of state vector errors to plot. Default is to use all which is likely what you want if the state space errors are standarized
#' @param t_index The bin indices to plot. Default is to use all bins
#' @param pch,ylab,xlab Arguments to override defaults set in the function
#' @param x_tick_loc,x_tick_mark \code{at} and \code{labels} arguments passed to \code{axis}
#' @param ... Arguments passed to plot
#'
#' @export
plot.fahrmeier_94_SpaceErrors = function(x, mod, cov_index = NA, t_index = NA,
                                         p_cex = par()$cex * .2, pch = 16,
                                         ylab = "Std. state space error",
                                         x_tick_loc = NA, x_tick_mark = NA,
                                         xlab = "Time", ...){
  bin_times = mod$times[-1]

  var_index = if(length(t_index) == 1 && is.na(cov_index)) seq_len(ncol(mod$state_vecs)) else cov_index
  res_std = x$residuals[, var_index, drop = F]
  n_vars = length(var_index)

  # assume equal distance
  delta_t = bin_times[2] - bin_times[1]
  delta_points = delta_t *.2 * (1:n_vars - n_vars/2)/n_vars

  # dummy plot
  use_custom_x_axis = !is.na(x_tick_loc) && ! is.na(x_tick_mark)
  t_index = if(length(t_index) == 1 && is.na(t_index)) seq_along(bin_times) else t_index
  plot(range(res_std) ~ c(min(bin_times[t_index]) + min(delta_points),
                          max(bin_times[t_index]) + max(delta_points)),
       type = "n", ylab = "Std. state space error", xlab = xlab,
       xaxt = ifelse(use_custom_x_axis, "n", "something"),
       ...)

  if(use_custom_x_axis)
    axis(1, x_tick_loc, x_tick_mark, lwd = par()$lwd)

  # add points
  for(i in t_index)
    points(bin_times[i] + delta_points, res_std[i, ], pch = pch, cex = p_cex)

  # add 95% conf
  abline(h = c(-1, 1) * 1.96, lty = 2)
}
