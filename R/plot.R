#' @title Plots for ddhazard Object
#' @description Plot of estimated state space variables from a \code{\link{ddhazard}} fit.
#'
#' @param x result of \code{\link{ddhazard}} call.
#' @param type type of plot. Currently, only \code{"cov"} is available for plot of the state space parameters.
#' @param plot_type the \code{type} argument passed to \code{plot}.
#' @param cov_index the index (indices) of the state space parameter(s) to plot.
#' @param add \code{FALSE} if you want to make a new plot.
#' @param xlab,ylab,ylim,col arguments to override defaults set in the function.
#' @param do_alter_mfcol \code{TRUE} if the function should alter \code{par(mfcol)} in case that \code{cov_index} has more than one element.
#' @param level level (fraction) for confidence bounds.
#' @param ddhazard_boot object from a \code{\link{ddhazard_boot}} call which confidence bounds will be based on and where bootstrap samples will be printed with a transparent color.
#' @param ... arguments passed to \code{\link{plot.default}} or \code{lines} depending on the value of \code{add}.
#'
#' @details
#' Creates a plot of state variables or adds state variables to a plot with indices \code{cov_index}. Pointwise 1.96 std. confidence intervals are provided with the smoothed co-variance matrices from the fit.
#'
#' @importFrom graphics matplot matpoints
#'
#' @examples
#'library(dynamichazard)
#'fit <- ddhazard(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
#'  Q_0 = diag(1, 2), Q = diag(1e-4, 2), by = 50,
#'  control = ddhazard_control(method = "GMA"))
#'plot(fit)
#'plot(fit, cov_index = 2)
#'
#' @importFrom graphics par plot lines
#' @importFrom grDevices col2rgb rgb
#' @export
plot.ddhazard = function(x, xlab = "Time",
                             ylab = "Hazard",
                             type = "cov", plot_type = "l", cov_index, ylim,
                             col = "black", add = F, do_alter_mfcol = T,
                             level = 0.95,
                             ddhazard_boot, ...){
  if(!x$model %in% c("logit", "cloglog", exp_model_names))
    stop("Functions for model '", x$model, "' is not implemented")

  missing_boot <- missing(ddhazard_boot)
  if(!missing_boot){
    # Removed failed estimates
    ddhazard_boot$t <- ddhazard_boot$t[!apply(ddhazard_boot$t, 1, function(x) any(is.na(x))), ]
  }

  if(type == "cov"){
    if(missing(cov_index)){
      n_cov <- dim(x$state_vecs)[2] / x$order
      if(n_cov > 0){
        cov_index <- 1:min(9, n_cov)
      } else
        stop("plot.ddhazard called with no time varying effects")
    }

    n_plots <- length(cov_index)
    if(!add && do_alter_mfcol && n_plots > 1){
      org_mfcol <- par()$mfcol
      on.exit(par(mfcol = org_mfcol))
      par(mfcol =
            if(n_plots <= 2) c(2,1) else
              if(n_plots <= 4) c(2,2) else
                if(n_plots <= 6) c(2,3) else
                  c(3,3))
    }

    for(i in cov_index){
      ylab_to_use <- if(missing(ylab)) colnames(x$state_vecs)[i] else ylab

      if(missing_boot){
        fac <- qnorm(.5 - level / 2)
        lb = x$state_vecs[, i] + fac * sqrt(x$state_vars[i, i, ])
        ub = x$state_vecs[, i] - fac * sqrt(x$state_vars[i, i, ])

      } else {
        boot_ests <- ddhazard_boot$t[, (i - 1) * nrow(x$state_vecs) + 1:nrow(x$state_vecs)]
        R <- nrow(boot_ests)
        rk_mat <- apply(boot_ests, 2, rank)

        # Find lower bound
        frac_n_weights <- get_frac_n_weights(R = R, a = .5 - level / 2)
        qs <- apply(rk_mat, 2, function(x){
          k <- frac_n_weights$k
          c(which(x == k), which(x == k + 1))})

        lb <- sapply(1:ncol(qs), function(x)
          frac_n_weights$w_k * boot_ests[qs[1, x], x] +
            frac_n_weights$w_k_p_1 * boot_ests[qs[2, x], x])

        # Do the same for the upper bound
        frac_n_weights <- get_frac_n_weights(R = R, a = .5 + level / 2)
        qs <- apply(rk_mat, 2, function(x){
          k <- frac_n_weights$k
          c(which(x == k), which(x == k + 1))})

        ub <- sapply(1:ncol(qs), function(x)
          frac_n_weights$w_k * boot_ests[qs[1, x], x] +
            frac_n_weights$w_k_p_1 * boot_ests[qs[2, x], x])
      }

      ylim_to_use <- if(missing(ylim)) range(lb, ub) else ylim

      if(!add){
        plot(x$times, x$state_vecs[, i], type = plot_type,
             ylim = ylim_to_use, xlab = xlab, ylab = ylab_to_use, col = col, ...)

      } else {
        lines(x$times, x$state_vecs[, i], col = col, ...)

      }

      lines(x$times, lb, lty = 2, col = col)
      lines(x$times, ub, lty = 2, col = col)

      if(!missing_boot){
        max_print = 50
        if(nrow(boot_ests) > max_print){
          if(i == cov_index[1])
            message("Only plotting ", max_print, " of the boot sample estimates")
          boot_ests <- boot_ests[sample.int(nrow(boot_ests), size = max_print, replace = F), ]
        }

        new_col <- c(col2rgb(col = col))
        new_col <- rgb(new_col[1], new_col[2], new_col[3], .05)
        for(i in 1:nrow(boot_ests)){
          lines(x$times, boot_ests[i, ], col = new_col)
        }
      }
    }

    return(invisible())
  }

  stop("Type '", type, "' is not implemented for plot.ddhazard")
}

#' @title State Space Error Plot
#' @description Plot function for state space errors from \code{\link{ddhazard}} fit.
#'
#' @param x result of \code{\link[=residuals.ddhazard]{residuals}} with a `type` argument which yields state space errors.
#' @param mod the \code{\link{ddhazard}} result used in the \code{\link[=residuals.ddhazard]{residuals}} call.
#' @param p_cex \code{cex} argument for the points
#' @param cov_index the indices of state vector errors to plot. Default is to use all.
#' @param t_index the bin indices to plot. Default is to use all bins.
#' @param pch,ylab,xlab arguments to override defaults set in the function.
#' @param x_tick_loc,x_tick_mark \code{at} and \code{labels} arguments passed to \code{axis}.
#' @param ... arguments passed to \code{\link{plot.default}}.
#'
#' @importFrom graphics abline axis par plot points
#' @export
plot.ddhazard_space_errors = function(x, mod, cov_index = NA, t_index = NA,
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
       type = "n", ylab = ylab, xlab = xlab,
       xaxt = ifelse(use_custom_x_axis, "n", "s"),
       ...)

  if(use_custom_x_axis)
    axis(1, x_tick_loc, x_tick_mark, lwd = par()$lwd)

  # add points
  for(i in t_index)
    points(bin_times[i] + delta_points, res_std[i, ], pch = pch, cex = p_cex)

  # add 95% conf
  abline(h = c(-1, 1) * 1.96, lty = 2)

  invisible()
}

#' @title Plot of Clouds from a PF_clouds Object
#' @description
#' Plots mean curve along with quantiles through time for the forward, backward or smoothed clouds.
#'
#' @param x an object of class \code{PF_clouds}.
#' @param y unused.
#' @param type parameter to specify which cloud to plot.
#' @param ylim \code{ylim} passed to \code{\link{matplot}}.
#' @param add \code{TRUE} if a new plot should not be made.
#' @param qlvls vector of quantile levels to be plotted.
#' @param pch \code{pch} argument for the quantile points.
#' @param lty \code{lty} argument for the mean curves.
#' @param col \code{col} argument to \code{\link{matplot}} and
#' \code{\link{matpoints}} or \code{\link{matlines}}.
#' @param ... unused.
#' @param cov_index indices of the state vector to plot. All are plotted if
#' this argument is omitted.
#' @param qtype character specifying how to show quantiles. Either
#' \code{"points"} for crosses or \code{"lines"} for dashed lines.
#'
#' @return
#' List with quantile levels and mean curve.
#'
#' @importFrom graphics matplot matpoints matlines
#' @export
plot.PF_clouds <- function(
  x, y,
  type = c("smoothed_clouds", "forward_clouds", "backward_clouds"),
  ylim, add = FALSE, qlvls = c(.05, .5, .95), pch = 4, lty = 1, col, ...,
  cov_index, qtype = c("points", "lines")){
  type <- type[1]
  these_clouds <- x[[type]]
  if(missing(cov_index))
    cov_index <- seq_len(dim(these_clouds[[1]]$states)[1])

  qtype <- qtype[1L]
  stopifnot(qtype %in% c("points", "lines"))

  #####
  # find means
  .mean <- get_cloud_means(x, cov_index = cov_index, type = type)

  #####
  # find quantiles
  if(length(qlvls) > 0){
    qs <- get_cloud_quantiles(
      x, cov_index = cov_index, type = type, qlvls = qlvls)

  } else
    qs <- NULL

  #####
  # plot
  .x <- 1:nrow(.mean) + if(type == "forward_clouds") -1 else 0
  if(missing(ylim))
    ylim <- range(qs, .mean, na.rm = TRUE) # can have NA if we dont have a
                                           # weighted value below or above qlvl

  if(missing(col))
    col <- 1:ncol(.mean)

  matplot(.x, .mean, ylim = ylim, type = "l", lty = lty, add = add, col = col,
          xlab = "Time", ylab = "State vector", ...)

  if(length(qs) > 0){
    for(i in 1:dim(qs)[2]){
      tmp <- qs[, i, ]
      if(is.null(dim(tmp)))
        tmp <- matrix(tmp, ncol = length(tmp))
      switch (qtype,
        points = matpoints(.x, t(tmp), pch = pch, col = col[i]),
        lines  = matlines(.x, t(tmp), lty = 3, col = col[i]),
        stop(sQuote(qtype), " not implemented"))
    }
  }

  invisible(list(mean = .mean, qs = qs))
}

#' @title Plot for a PF_EM Object
#' @description
#' Short hand to call \code{\link{plot.PF_clouds}}.
#'
#' @param x an object of class \code{PF_EM}.
#' @param y unused.
#' @param ... arguments to \code{\link{plot.PF_clouds}}.
#'
#' @return
#' See \code{\link{plot.PF_clouds}}
#'
#' @export
plot.PF_EM <- function(x, y, ...){
  cl <- match.call()
  cl[[1]] <- quote(plot)
  cl$x <- bquote(.(substitute(x))$clouds)
  invisible(eval(cl, parent.frame()))
}

#' @describeIn  ddsurvcurve  method for plotting survival curve.
#' @section plot.ddsurvcurve:
#' Returns the same as \code{lines.ddsurvcurve}.
#' @export
plot.ddsurvcurve <- function(x, y, xlab = "Time", ylab = "Survival", ylim,
                             xaxs = "i", yaxs = "i", ...){
  if(missing(ylim))
    ylim <- c(0, 1.04)
  if(yaxs == "i" && ylim[2] == 1)
    ylim[2] <- 1 + .04 * diff(ylim)

  plot(range(x$start, x$time), c(0, 1), type = "n", xlab = xlab, ylab = ylab,
       ylim = ylim, xaxs = xaxs, yaxs = yaxs)
  cl <- match.call()
  cl[[1L]] <- quote(lines)
  cl[c("xlab", "Survial")] <- NULL
  invisible(eval(cl, parent.frame()))
}

#' @describeIn ddsurvcurve Method for adding survival curve to a plot.
#' @section lines.ddsurvcurve:
#' Either returns the objects used in the call to \code{\link{segments}} for discrete
#' time hazard models, or the time points and survival function used to draw
#' the survival curve.
#' @importFrom graphics segments par
#' @export
lines.ddsurvcurve <- function(x, col = "Black", lty = 1, lwd = par()$lwd, ...){
  tstart <- c(x$start, x$time[-length(x$time)])
  tstop <- x$time

  if(x$family %in% c("logistic", "cloglog")){
    yold <- c(1, x$psurv[-length(x$psurv)])
    segments(tstart, yold, tstop, yold, col = col, lty = lty, lwd = lwd)
    segments(tstop, yold, tstop, x$psurv, col = col, lty = lty, lwd = lwd)

    return(invisible(list(tstart = tstart, yold = yold, psurv = x$psurv)))

  } else if(x$family == "exponential"){
    # find rates
    tdiff <- tstop  - tstart
    lambdas <- -log(1 - x$dhazard) / tdiff
    survprev <- c(1, x$psurv[-length(x$psurv)])
    bin_borders <- c(x$start, x$time)
    bin_borders[length(bin_borders)] <- bin_borders[length(bin_borders)] + 1e-8

    # assign survival function
    S <- function(z){
      idx <- findInterval(z, bin_borders, left.open = FALSE)
      if(any(idx <= 0, idx >= max(bin_borders)))
        stop("data points outside time range")

      f1 <- survprev[idx]
      f2 <- exp(-lambdas[idx] * (z - tstart[idx]))
      f1 * f2
    }

    z <- seq(tstart[1], tstop[length(tstop)], length.out = 1000)
    lines(z, S(z), col = col, lty = lty, lwd = lwd)

    return(invisible(list(z = z, S = S)))
  }

  stop(sQuote("lines.ddsurvcurve"), " is not implemented for family ",
       sQuote(x$family))
}
