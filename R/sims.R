# Simple function to make sure that we do not call rexp to many times when
# simulating
get_exp_draw <- function(cmpfun = TRUE) {
  out <- with(new.env(parent = .GlobalEnv), {
    is_first <- TRUE
    n_max <- n_draws <- 1e5
    n_cur <- NULL
    draws <- NULL

    function(n, re_draw = F){
      if(re_draw || is_first){
        n_max <<- n_draws
        n_cur <<- 1
        draws <<- rexp(n = n_draws, rate = 1)
        is_first <<- FALSE

        if(re_draw)
          return(NULL)
      }

      if(n == 1 && n_max > n_cur){
        n_cur <<- n_cur + 1
        return(draws[n_cur - 1])
      }

      # draw if needed
      if(n > n_max - n_cur){ # we forget about the - 1
        n_max <<- n + n_draws
        n_cur <<- 1
        draws <<- rexp(n = n + n_draws, rate = 1)
      }

      tmp <- n_cur
      n_cur <<- n_cur + n
      draws[tmp:(tmp + (n - 1))]
    }
  })

  if(cmpfun)
    compiler::cmpfun(out, options = list(optimize = 2, suppressAll = T)) else
      out
}

get_norm_draw <- function(cmpfun = TRUE) {
  out <- with(new.env(parent = .GlobalEnv), {
    is_first <- TRUE
    n_max <- n_draws <- 1e5
    n_cur <- NULL
    draws <- NULL

    function(n, re_draw = F){
      if(re_draw || is_first){
        n_max <<- n_draws
        n_cur <<- 1
        draws <<- rnorm(n = n_draws)
        is_first <<- FALSE

        if(re_draw)
          return(NULL)
      }

      if(n == 1 && n_max > n_cur){
        n_cur <<- n_cur + 1
        return(draws[n_cur - 1])
      }

      # Draw if needed
      if(n > n_max - n_cur){ # we forget about the - 1
        n_max <<- n + n_draws
        n_cur <<- 1
        draws <<- rnorm(n = n + n_draws)
      }

      tmp <- n_cur
      n_cur <<- n_cur + n
      draws[tmp:(tmp + (n - 1))]
    }
  })

  if(cmpfun)
    compiler::cmpfun(out, options = list(optimize = 2, suppressAll = T)) else
      out
}

get_unif_draw <- function(cmpfun = TRUE) {
  out <- with(new.env(parent = .GlobalEnv), {
    is_first <- TRUE
    n_max <- n_draws <- 1e5
    n_cur <- NULL
    draws <- NULL

    function(n, re_draw = F){
      if(re_draw || is_first){
        n_max <<- n_draws
        n_cur <<- 1
        draws <<- runif(n = n_draws)
        is_first <<- FALSE

        if(re_draw)
          return(NULL)
      }

      if(n == 1 && n_max > n_cur){
        n_cur <<- n_cur + 1
        return(draws[n_cur - 1])
      }

      # Draw if needed
      if(n > n_max - n_cur){ # we forget about the - 1
        n_max <<- n + n_draws
        n_cur <<- 1
        draws <<- runif(n = n + n_draws)
      }

      tmp <- n_cur
      n_cur <<- n_cur + n
      draws[tmp:(tmp + (n - 1))]
    }
  })

  if(cmpfun)
    compiler::cmpfun(out, options = list(optimize = 2, suppressAll = T)) else
      out
}

# n <- 5e4
# microbenchmark(
#   { f <- get_exp_draw(); for(i in 1:n) { f(1) }},
#   for(i in 1:n) { rexp(1,1) })

# define functions to simulate outcomes
test_sim_func_logit <- function(
  n_series, n_vars = 10L, t_0 = 0L, t_max = 10L, x_range = .1, x_mean = -.1,
  re_draw = T, beta_start = 3, intercept_start,
  sds = rep(1, n_vars + !missing(intercept_start)),
  is_fixed = c(), lambda = 1,
  tstart_sampl_func = function(t_0 = t_0, t_max = t_max) t_0, betas){
  cl <- match.call()
  cl[["linkfunc"]] <- "logit"
  cl[[1L]] <- bquote(
    environment(.(cl[[1L]]))$test_sim_func_discrete)
  eval(cl, parent.frame())
}

test_sim_func_discrete <- function(
  n_series, n_vars = 10L, t_0 = 0L, t_max = 10L, x_range = .1, x_mean = -.1,
  re_draw = T, beta_start = 3, intercept_start,
  sds = rep(1, n_vars + !missing(intercept_start)),
  is_fixed = c(), lambda = 1,
  tstart_sampl_func = function(t_0 = t_0, t_max = t_max) t_0, betas,
  linkfunc){
  # make output matrix
  n_row_max <- n_row_inc <- 10^5
  res <- matrix(NA_real_, nrow = n_row_inc, ncol = 4 + n_vars,
                dimnames = list(NULL, c("id", "tstart", "tstop", "event", paste0("x", 1:n_vars))))
  cur_row <- 1

  get_unif_draw <- get_unif_draw()
  get_exp_draw <- get_exp_draw()
  get_norm_draw <- get_norm_draw()

  if(re_draw){
    get_unif_draw(re_draw = T)
    get_exp_draw(re_draw = T)
    get_norm_draw(re_draw = T)
  }

  if(length(beta_start) == 1)
    beta_start <- rep(beta_start, n_vars)

  # draw betas
  if(missing(betas)){
    use_intercept <- !missing(intercept_start)
    betas <- matrix(get_norm_draw((t_max - t_0 + 1) * (n_vars + use_intercept)),
                    ncol = n_vars + use_intercept, nrow = t_max - t_0 + 1)
    betas <- t(t(betas) * sds)
    betas[1, ] <- if(use_intercept) c(intercept_start, beta_start) else beta_start
    betas <- apply(betas, 2, cumsum)

    betas[, is_fixed] <- matrix(rep(betas[1, is_fixed], nrow(betas)), byrow = T,
                                nrow = nrow(betas))
  } else
    use_intercept = ncol(betas) >= n_vars

  ceiler <- function(x)
    ceiling(x * 100) / 100
  x_adj <- - x_range / 2 + x_mean

  linkfunc <- switch(
    linkfunc,
    logit = function(x) 1 / (1 + exp(-x)),
    cloglog = function(x) -expm1(-exp(x)),
    stop(sQuote("linkfunc"), " not implemented"))

  # simulate
  for(id in 1:n_series){
    tstart <- tstop <- ceiler(tstart_sampl_func(t_0, t_max))
    interval_start <- ceiling(tstart)
    repeat{
      tstop <- ceiler(tstart + 1 / lambda * get_exp_draw(1) + 1)
      if(tstop >= t_max)
        tstop <- t_max

      x_vars <- x_range * get_unif_draw(n_vars) + x_adj
      l_x_vars <- if(use_intercept) c(1, x_vars) else x_vars

      tmp_t <- tstart
      while(tmp_t <= interval_start &&  interval_start < tstop) {
        event <- linkfunc((betas[interval_start + 2, ] %*% l_x_vars)[1, 1]) >
          get_unif_draw(1)

        interval_start <- interval_start + 1L
        if(event || interval_start >= t_max){
          tstop <- interval_start
          break
        }

        tmp_t <- tmp_t + 1
      }

      res[cur_row, ] <- c(id, tstart, tstop, event, x_vars)

      if(cur_row == n_row_max){
        n_row_max <- n_row_max + n_row_inc
        res = rbind(res, matrix(NA_real_, nrow = n_row_inc, ncol = 4 + n_vars))
      }
      cur_row <- cur_row + 1

      if(event || interval_start >= t_max)
        break

      tstart <- tstop
    }
  }

  list(res = as.data.frame(res[1:(cur_row - 1), ]), betas = betas)
}

# tmp_file <- tempfile()
# Rprof(tmp_file)
# tmp <- test_sim_func_logit(10^5)
# Rprof(NULL)
# summaryRprof(tmp_file)
# #
# sum(tmp$res[, "event"])

test_sim_func_exp <- function(
  n_series, n_vars = 10, t_0 = 0, t_max = 10, x_range = 1, x_mean = 0,
  re_draw = T, beta_start = 1, intercept_start,
  sds = rep(1, n_vars + !missing(intercept_start)),
  is_fixed = c(), lambda = 1,
  tstart_sampl_func = function(t_0 = t_0, t_max = t_max)
    t_0){
  # make output matrix
  n_row_max <- n_row_inc <- 10^5
  res <- matrix(NA_real_, nrow = n_row_inc, ncol = 4 + n_vars,
                dimnames = list(NULL, c("id", "tstart", "tstop", "event",
                                        paste0("x", 1:n_vars))))
  cur_row <- 1

  get_unif_draw <- get_unif_draw()
  get_exp_draw <- get_exp_draw()
  get_norm_draw <- get_norm_draw()

  if(re_draw){
    get_unif_draw(re_draw = T)
    get_exp_draw(re_draw = T)
    get_norm_draw(re_draw = T)
  }

  if(length(beta_start) == 1)
    beta_start <- rep(beta_start, n_vars)

  # draw betas
  use_intercept <- !missing(intercept_start)
  betas <- matrix(get_norm_draw((t_max - t_0 + 1) * (n_vars + use_intercept)),
                  ncol = n_vars + use_intercept, nrow = t_max - t_0 + 1)
  betas <- t(t(betas) * sds)
  betas[1, ] <-
    if(use_intercept) c(intercept_start, beta_start) else beta_start
  betas <- apply(betas, 2, cumsum)

  betas[, is_fixed] <- matrix(rep(betas[1, is_fixed], nrow(betas)), byrow = T,
                              nrow = nrow(betas))

  ceiler <- function(x, level=1) round(x + 5*10^(-level-1), level)

  # simulate
  mean_term <- - x_range / 2 + x_mean

  for(id in 1:n_series){
    tstart <- tstop <-  ceiler(tstart_sampl_func(t_0, t_max), 2)
    repeat{
      tstop <- ceiler(tstart + get_exp_draw(1)/lambda, 2)

      x_vars <- x_range * get_unif_draw(n_vars) + mean_term
      l_x_vars <- if(use_intercept) c(1, x_vars) else x_vars

      tmp_t <- tstart
      while(tmp_t < tstop && tmp_t < t_max){
        delta_max <- min(ceiling(tmp_t + 1e-14), tstop) - tmp_t
        new_time <-
          get_exp_draw(1) /
          exp(drop(betas[floor(tmp_t - t_0) + 2, ] %*% l_x_vars))
        event <- new_time <= delta_max
        if(event){
          tstop <- ceiler(new_time + tmp_t, 14)
          break
        }

        tmp_t <- ceiling(tmp_t + 1e-14)
      }

      res[cur_row, ] <- c(id, tstart, tstop, event, x_vars)

      if(cur_row == n_row_max){
        n_row_max <- n_row_max + n_row_inc
        res = rbind(res, matrix(NA_real_, nrow = n_row_inc, ncol = 4 + n_vars))
      }
      cur_row <- cur_row + 1

      if(event || tstop >= t_max)
        break

      tstart <- tstop
    }
  }

  list(res = as.data.frame(res[1:(cur_row - 1), ]), betas = betas)
}
