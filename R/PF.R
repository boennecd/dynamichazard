#' @export
logLik.PF_clouds <- function(object){
  sum(tail(
    sapply(lapply(object$forward_clouds,
                  "[[", "log_unnormalized_weights"),
           function(x){
             .max <- max(x)
             log(sum(exp(x - .max))) + .max - log(length(x))
          }),
    -1))
}

PF_effective_sample_size <- function(object){
  sapply(object, function(x){
    sapply(lapply(x, "[[", "weights"), function(z) 1 / sum(z^2))
  })
}

# TODO: test
# TODO: remove export
#' @export
PF_EM <- function(
  formula, data,
  model = "logit",
  by, max_T, id,
  a_0, Q_0, Q = Q_0,
  order = 1,
  control = list(),
  verbose = F
){
  #####
  # Checks
  if(order != 1) # TODO: test
    stop(sQuote('order'), " not equal to 1 is not supported")

  if(!model %in% "logit") # TODO: test
    stop(sQuote('model'), " is not supported")

  if(missing(id)){
    if(verbose)
      warning("You did not parse and Id argument")
    id = 1:nrow(data)
  }

  #####
  # FInd design matrix

  X_Y = get_design_matrix(formula, data)
  n_params = ncol(X_Y$X)

  if(length(X_Y$fixed_terms) > 0) # TODO: test
    stop("Fixed terms are not supported")


  #####
  # Set control variables
  control_default <- list(
    eps = 1e-2, forward_backward_ESS_threshold = NULL,
    method = "AUX_normal_approx_w_particles",
    trace = 0, n_max = 25,
    n_threads = getOption("ddhazard_max_threads"))

  if(any(is.na(control_match <- match(names(control), names(control_default)))))
    stop("These control parameters are not recognized: ",
         paste0(names(control)[is.na(control_match)], collapse = "\t"))

  control_default[control_match] <- control
  control <- control_default

}


.PF_EM <- function(
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
      plot(clouds, type = "forward_clouds", add = TRUE, qlvls = c(), lty = 2)
      plot(clouds, type = "backward_clouds", add = TRUE, qlvls = c(), lty = 3)

      cat("Effective sample sizes are:\n")
      print(effective_sample_size <- PF_effective_sample_size(clouds))
    }

    #####
    # Update parameters
    sum_stats <- compute_summary_stats(clouds)
    a_0 <- drop(sum_stats[[1]]$E_xs)
    Q <- matrix(0., length(a_0), length(a_0))
    for(j in 2:length(sum_stats))
      Q <- Q + sum_stats[[j]]$E_x_less_x_less_one_outers
    Q <- Q / (length(sum_stats) - 1)

    args$a_0 <- a_0
    args$Q <- Q

    #####
    # Compute log likelihood and check for convergernce
    log_like <- logLik(clouds)
    log_likes[i] <- log_like

    if(trace > 0)
      cat("The log likelihood in iteration ", i, " is ", log_like,
          ". Largest log likelihood before this iteration is ", log_like_max, "\n", sep = "")

    if(log_like < log_like_max)
      warning("Likelihood decreased in iteration ", i, " compared to maximum likelihood seen so far.")

    if(has_converged <- log_like_max < log_like && log_like - log_like_max < eps)
      break
  }

  if(!has_converged)
    warning("Method did not converge.")

  if(!exists("effective_sample_size", envir = environment()))
    effective_sample_size <- PF_effective_sample_size(clouds)

  return(structure(list(
    call = cl,
    clouds = clouds,
    a_0 = a_0,
    Q = Q,
    F = args$F,
    summary_stats = sum_stats,
    log_likes = log_likes[1:i],
    n_iter = i,
    effective_sample_size = effective_sample_size,
    seed = seed),
    class = "PF_EM"))
}

# TODO: test
# TODO: move to plot.R
#' @export
plot.PF_clouds <- function(
  x, y,
  type = c("smoothed_clouds", "forward_clouds", "backward_clouds"),
  ylim, add = FALSE, qlvls = c(.05, .5, .95), pch = 4, lty = 1, ..., cov_index){
  type <- type[1]
  these_clouds <- x[[type]]
  if(missing(cov_index))
    cov_index <- seq_len(dim(these_clouds[[1]]$states)[1])

  #####
  # Find means
  .mean <- t(sapply(these_clouds, function(row){
    colSums(t(row$states[cov_index, , drop = FALSE]) * drop(row$weights))
  }))

  #####
  # Find quantiles
  if(length(qlvls) > 0){
    qs <- lapply(these_clouds, function(row){
      out <- apply(row$states[cov_index, , drop = FALSE], 1, function(x){
        ord <- order(x)
        wg_cumsum <- cumsum(row$weights[ord])
        idx <- ord[sapply(qlvls, function(q) max(which(wg_cumsum < q)))]
        x[idx]
      })

      if(is.null(dim(out)))
        out <- matrix(out, ncol = length(out))

      out
      })
    qs <- simplify2array(qs)
  } else
    qs <- NULL

  #####
  # Plot
  .x <- 1:nrow(.mean) + if(type == "forward_clouds") -1 else 0
  if(missing(ylim))
    ylim <- range(qs, .mean, na.rm = TRUE) # can have NA if we dont have a
                                           # weighted value below or above
                                           # qlvl
  matplot(.x, .mean, ylim = ylim, type = "l", lty = lty, add = add,
          xlab = "Time", ylab = "State vector", ...)
  if(length(qs) > 0){
    for(i in 1:dim(qs)[2]){
      tmp <- qs[, i, ]
      if(is.null(dim(tmp)))
        tmp <- matrix(tmp, ncol = length(tmp))
      matpoints(.x, t(tmp), pch = pch, col = i)
    }
  }

  invisible(list(mean = .mean, qs = qs))
}

# TODO: test
# TODO: move to plot.R
#' @export
plot.PF_EM <- function(x, y, ...){
  invisible(do.call(plot, c(list(x = x$clouds), list(...))))
}

