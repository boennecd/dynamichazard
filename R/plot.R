#' Plots to illustrate the estimate state space variables from a ddhazardfit
#' @export
plot.fahrmeier_94 = function(object, new_data, xlab = "Time",
                             ylab = "Hazard", method = "delta",
                             type, plot_type = "l", cov_index, ylim, ...){
  if(!object$model %in% c("logit"))
    stop("Functions for model '", object$model, "' is not implemented")

  if(!missing(type) && type == "cov"){
    if(missing(cov_index))
      stop("Need cov index when type is equal to ", type)

    if(missing(ylab))
      ylab = colnames(object$a_t_d_s)[cov_index]

    lb = object$a_t_d_s[, cov_index] - 1.96 * sqrt(object$V_t_d_s[cov_index, cov_index, ])
    ub = object$a_t_d_s[, cov_index] + 1.96 * sqrt(object$V_t_d_s[cov_index, cov_index, ])

    if(missing(ylim))
      ylim = range(lb, ub)

    plot(object$times, object$a_t_d_s[, cov_index], type = plot_type,
         ylim = ylim, xlab = xlab, ylab = ylab, ...)
    lines(object$times, lb, lty = 2)
    lines(object$times, ub, lty = 2)

    return(invisible())
  }

  predict_ = predict.fahrmeier_94(object, new_data, method = method)
  hazard = predict_$fits / c(1, diff(object$times))

  plot(predict_$times, hazard, xlab = xlab, ylab = ylab, type = plot_type, ...)
  lines(predict_$times, predict_$lbs, lty = 2)
  lines(predict_$times, predict_$ubs, lty = 2)
}

#' Plot function for standadized statespace errors from ddhazardfit
#' @export
plot.fahrmeier_94_stdSpaceErrors = function(object, mod, vars_ = NA, t_index = NA,
                                            p_cex = par()$cex * .2, pch = 16,
                                            ylab = "Std. state space error",
                                            x_tick_loc = NA, x_tick_mark = NA,
                                            xlab = "Time", ...){
  bin_times = mod$times[-1]

  var_index = if(length(t_index) == 1 && is.na(vars_)) seq_len(ncol(mod$a_t_d_s)) else vars_
  res_std = object[, var_index, drop = F]
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
