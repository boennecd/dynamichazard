logLik.PF_clouds <- function(object){
  sum(tail(
    sapply(lapply(object$forward_clouds,
                  "[[", "log_unnormalized_weights"), mean),
    -1))
}

PF_effective_sample_size <- function(object){
  sapply(object, function(x){
    sapply(lapply(x, "[[", "weights"), function(z) 1 / sum(z^2))
  })
}

PF_EM <- function(
  n_fixed_terms_in_state_vec,
  X,
  fixed_terms,
  tstart,
  tstop,
  Q_0,
  Q,
  a_0 ,
  risk_obj,
  n_max,
  order,
  n_threads,
  N_fw_n_bw,
  N_smooth,
  N_first,
  eps,
  forward_backward_ESS_threshold = NULL,
  trace = 0,
  method = "AUX_normal_approx_w_particles",
  seed){
  if(order != 1) # TODO: test
    stop(sQuote('order'), " not equal to 1 is not supported")

  if(length(X_Y$fixed_terms) > 0) # TODO: test
    stop("Fixed terms are not supported")

  cl <- match.call()
  args <- as.list(cl)
  n_vars <- nrow(X)
  args[[1]] <- NULL
  args[["Q_tilde"]] <- diag(0, n_vars)
  args[["F"]] <- diag(1, n_vars)
  args[["debug"]] <- max(0, trace - 1)
  args[c("trace", "eps", "seed")] <- NULL

  if(missing(seed)){
    seed <- .Random.seed
  }

  log_likes <- rep(NA_real_, n_max)
  log_like <- log_like_max <- -Inf
  for(i in 1:n_max){
    if(trace > 0){
      if(i != 1)
        cat("\n\n\n")
      cat("#######################\nStarting EM iteration", i, "\n")

      cat("a_0 is:\n")
      print(eval(args$a_0, environment()))

      cat("chol(Q) is:\n")
      print(chol(eval(args$Q, environment())))
    }
    log_like_old <- log_like
    log_like_max <- max(log_like, log_like_max)

    #####
    # Find clouds
    set.seed(seed)
    clouds <- do.call(PF_smooth, args)

    if(trace > 0){
      cat("Plotting state vector mean and quantiles for iteration", i, "\n")
      plot(clouds, main = paste0("EM iteration ", i))
    }

    #####
    # Update parameters
    sum_stats <- compute_summary_stats(clouds)
    a_0 <- drop(sum_stats[[1]]$E_xs)
    Q <- matrix(0., length(a_0), length(a_0))
    for(j in 1:length(sum_stats))
      Q <- Q + sum_stats[[j]]$E_x_less_x_less_one_outers
    Q <- Q / length(sum_stats)

    args$a_0 <- a_0
    args$Q <- Q

    #####
    # Compute log likelihood and check for convergernce
    log_like <- logLik(clouds)
    log_likes[i] <- log_like

    if(trace > 0)
      cat("The log likelihood in iteration", i, "is", log_like, "\n")

    if(log_like < log_like_max)
      warning("Likelihood decreased in iteration ", i, " compared to maximum likelihood seen so far.")

    if(has_converged <- log_like_max < log_like && log_like - log_like_max < eps)
      break
  }

  if(!has_converged)
    warning("Method did not converge.")

  return(structure(list(
    clouds = clouds,
    a_0 = a_0,
    Q = Q,
    F = args$F,
    summary_stats = sum_stats,
    log_likes = log_likes[1:i],
    n_iter = i,
    effective_sample_size = PF_effective_sample_size(clouds),
    seed = seed),
    class = "PF_EM"))
}

# TODO: test
plot.PF_clouds <- function(
  x, y,
  type = c("smoothed_clouds", "forward_clouds", "backward_clouds"),
  ylim, add = FALSE, qlvls = c(.025, .5, .975), ...){
  type <- type[1]
  these_clouds <- x[[type]]

  #####
  # Find means
  .mean <- t(sapply(these_clouds, function(row){
    colSums(t(row$states) * drop(row$weights))
  }))

  #####
  # Find quantiles
  qs <- lapply(these_clouds, function(row){
    apply(row$states, 1, function(x){
      ord <- order(x)
      wg_cumsum <- cumsum(row$weights[ord])
      idx <- ord[sapply(qlvls, function(q) max(which(wg_cumsum < q)))]
      x[idx]
    })
  })
  qs <- simplify2array(qs)

  #####
  # Plot
  .x <- 1:nrow(.mean) + if(type == "forward_clouds") -1 else 0
  if(missing(ylim))
    ylim <- range(qs, .mean, na.rm = TRUE) # can have NA if we dont have a
                                           # weighted value below or above
                                           # qlvl
  matplot(.x, .mean, ylim = ylim, type = "l", lty = 1, add = add,
          xlab = "Time", ylab = "State vector", ...)
  for(i in 1:dim(qs)[2]){
    matpoints(.x, t(qs[, i, ]), pch = 4, col = i)
  }

  invisible(list(mean = .mean, qs = qs))
}

# TODO: test
plot.PF_EM <- function(x, y, ...){
  invisible(do.call(plot, c(list(x = x$clouds), list(...))))
}

