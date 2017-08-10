# function to test that expression does not throw an error
expect_no_error = function(expr, env = parent.frame()){
  call <- match.call()
  eval(bquote(expect_error(.(call[[2]]), NA)), envir = env)
}

# expect_no_error(1 / "a")
# expect_no_error(1 / 1)

# Load data example to play arround with and validate against
get_head_neck_cancer_data <- function(){
  is_censored = c(6, 27, 34, 36, 42, 46, 48:51,
                  51 + c(15, 30:28, 33, 35:37, 39, 40, 42:45))
  head_neck_cancer = data.frame(
    id = 1:96,
    start = rep(0, 96),
    stop = c(
      1, 2, 2, rep(3, 6), 4, 4, rep(5, 8),
      rep(6, 7), 7, 8, 8, 8,
      9, 9, 10, 10, 10, 11, 14, 14, 14, 15, 18, 18, 20, 20, 37,
      37, 38, 41, 45, 47, 47,

      2, 2, 3, rep(4, 4), rep(5, 5), rep(6, 5),
      7, 7, 7, 9, 10, 11, 12, 15, 16, 18, 18, 18, 21,
      21, 24, 25, 27, 36, 41, 44, 52, 54, 59, 59, 63, 67, 71, 76),
    event = !(1:96 %in% is_censored),
    group = factor(c(rep(1, 45 + 6), rep(2, 45))))

  head_neck_cancer$group = factor(head_neck_cancer$group, levels = c(2, 1))

  head_neck_cancer
}

# Simple function to make sure that we do not call rexp to many times
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

      # Draw if needed
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
#   { get_exp_draw()(n); for(i in 1:n) { }},
#   for(i in 1:n) { rexp(1,1) })

# Define functions to simulate outcomes
test_sim_func_logit <- function(n_series, n_vars = 10L, t_0 = 0L, t_max = 10L, x_range = .1, x_mean = -.1,
                                re_draw = T, beta_start = 3, intercept_start,
                                sds = rep(1, n_vars + !missing(intercept_start)),
                                is_fixed = c(), lambda = 1,
                                tstart_sampl_func = function(t_0 = t_0, t_max = t_max)
                                  t_0,
                                betas){
  # Make output matrix
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
  } else {
    use_intercept = ncol(betas) >= n_vars
  }

  ceiler <- function(x) ceiling(x * 100) / 100
  x_adj <- - x_range / 2 + x_mean

  # Simulate
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
        exp_eta <- exp((betas[interval_start + 2, ] %*% l_x_vars)[1, 1])
        event <- exp_eta / (1 + exp_eta) > get_unif_draw(1)

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

test_sim_func_exp <- function(n_series, n_vars = 10, t_0 = 0, t_max = 10, x_range = 1, x_mean = 0,
                              re_draw = T, beta_start = 1, intercept_start,
                              sds = rep(1, n_vars + !missing(intercept_start)),
                              is_fixed = c(), lambda = 1,
                              tstart_sampl_func = function(t_0 = t_0, t_max = t_max)
                                t_0){
  # Make output matrix
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
  use_intercept <- !missing(intercept_start)
  betas <- matrix(get_norm_draw((t_max - t_0 + 1) * (n_vars + use_intercept)),
                  ncol = n_vars + use_intercept, nrow = t_max - t_0 + 1)
  betas <- t(t(betas) * sds)
  betas[1, ] <- if(use_intercept) c(intercept_start, beta_start) else beta_start
  betas <- apply(betas, 2, cumsum)

  betas[, is_fixed] <- matrix(rep(betas[1, is_fixed], nrow(betas)), byrow = T,
                              nrow = nrow(betas))

  ceiler <- function(x, level=1) round(x + 5*10^(-level-1), level)

  # Simulate
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
        new_time <- get_exp_draw(1) / exp(drop(betas[floor(tmp_t - t_0) + 2, ] %*% l_x_vars))
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

########
# PBC data set from survival with the timevariying covariates
# See: https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf

# Remove a globally defined pbc dataset if it is present
get_pbc2_data <- function(){
  pbc <- survival::pbc
  pbcseq <- survival::pbcseq

  # Avoid notes with CRAN tests
  id <- NULL
  sex <- NULL
  time <- NULL
  status <- NULL
  edema <- NULL
  age <- NULL
  event <- NULL
  tdc <- NULL
  day <- NULL
  albumin <- NULL
  protime <- NULL
  bili <- NULL

  temp <- subset(pbc, id <= 312, select=c(id, sex, time, status, edema, age))
  pbc2 <- survival::tmerge(
    temp, temp, id=id, death = event(time, status))
  pbc2 <- survival::tmerge(
    pbc2, pbcseq, id=id, albumin = tdc(day, albumin),
    protime = tdc(day, protime), bili = tdc(day, bili))
  pbc2 <- pbc2[, c("id", "tstart", "tstop", "death", "sex", "edema",
                   "age", "albumin", "protime", "bili")]

  pbc2
}

