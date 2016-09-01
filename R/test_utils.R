# function to test that expression does not throw an error
expect_no_error = function(expr, env = parent.frame()){
  call <- match.call()
  eval(bquote(expect_error(.(call[[2]]), NA)), envir = env)
}

# expect_no_error(1 / "a")
# expect_no_error(1 / 1)

# Load data example to play arround with and validate against
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

rm(is_censored)

head_neck_cancer$group = factor(head_neck_cancer$group, levels = c(2, 1))

# Simple function to make sure that we do not call rexp to many times
get_exp_draw <- with(new.env(), {
  n_max <- n_draws <- 10^5
  n_cur <- 1
  draws <- rexp(n = n_draws, rate = 1)
  #env_p <- environment()

  function(n, re_draw = F){
    if(re_draw){
      n_max <<- n_draws
      n_cur <<- 1
      draws <<- rexp(n = n_draws, rate = 1)

      return(NULL)
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

get_norm_draw <- with(new.env(), {
  n_max <- n_draws <- 10^5
  n_cur <- 1
  draws <- rnorm(n = n_draws)
  #env_p <- environment()

  function(n, re_draw = F){
    if(re_draw){
      n_max <<- n_draws
      n_cur <<- 1
      draws <<- rnorm(n = n_draws)

      return(NULL)
    }

    # Draw if needed
    if(n > n_max - n_cur){ # we forget about the - 1
      n_max <<- n + n_draws
      n_cur <<- 1
      draws <<- rnorm(n = n + n_draws, rate = 1)
    }

    tmp <- n_cur
    n_cur <<- n_cur + n
    draws[tmp:(tmp + (n - 1))]
  }
})
get_unif_draw <- with(new.env(), {
  n_max <- n_draws <- 10^5
  n_cur <- 1
  draws <- runif(n_draws)
  #env_p <- environment()

  function(n, re_draw = F){
    if(re_draw){
      n_max <<- n_draws
      n_cur <<- 1
      draws <<- runif(n = n_draws)

      return(NULL)
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

# library(microbenchmark)
# microbenchmark(get_exp_draw(1), rexp(1,1))
# microbenchmark(get_exp_draw(010), rexp(100,1))

get_exp_draw = compiler::cmpfun(get_exp_draw, options = list(
  optimize = 3, suppressAll = T))
get_unif_draw = compiler::cmpfun(get_unif_draw, options = list(
  optimize = 3, suppressAll = T))
get_norm_draw = compiler::cmpfun(get_norm_draw, options = list(
  optimize = 3, suppressAll = T))

# microbenchmark(get_exp_draw(1), rexp(1,1))
# microbenchmark(get_exp_draw(100), rexp(100,1))
# microbenchmark(get_unif_draw(1), runif(1))
# microbenchmark(get_unif_draw(100), runif(100))

# Define function to simulate outcomes
test_sim_func_logit <- function(n_series, n_vars = 10, t_0 = 0, t_max = 10, x_range = .1, x_mean = -.1,
                                re_draw = T, beta_start = 3){
  # Make output matrix
  n_row_max <- n_row_inc <- 10^5
  res <- matrix(NA_real_, nrow = n_row_inc, ncol = 4 + n_vars,
                dimnames = list(NULL, c("id", "tstart", "tstop", "event", paste0("x", 1:n_vars))))
  cur_row <- 1

  if(re_draw){
    get_unif_draw(re_draw = T)
    get_exp_draw(re_draw = T)
    get_norm_draw(re_draw = T)
  }

  # draw betas
  betas <- matrix(get_norm_draw((t_max - t_0) * n_vars), ncol = n_vars, nrow = t_max - t_0)
  betas[1, ] <- beta_start
  betas <- apply(betas, 2, cumsum)

  # Simulate
  for(id in 1:n_series){
    tstart <- tstop <- t_0
    repeat{
      tstop <- tstart + get_exp_draw(1) + 1

      x_vars <- x_range * get_unif_draw(n_vars) - x_range / 2 + x_mean
      tmp_t <- tstart
      while(tmp_t < tstop && tmp_t < t_max){
        exp_eta <- exp((betas[floor(tmp_t) + 1 + t_0, ] %*% x_vars)[1, 1])
        event <- exp_eta / (1 + exp_eta) > get_unif_draw(1)
        if(event){
          tstop <- tmp_t + 1
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

      if(event || tstop > t_max)
        break

      tstart <- tstop
    }
  }

  list(res = res[1:(cur_row - 1), ], betas = betas)
}

test_sim_func_logit = compiler::cmpfun(test_sim_func_logit)

# tmp_file <- tempfile()
# Rprof(tmp_file)
# tmp <- test_sim_func_logit(10^5)
# Rprof(NULL)
# summaryRprof(tmp_file)
# #
# sum(tmp$res[, "event"])

# hint use this regexp '\r\n\s*\[\d+\]\s' and replace with '\r\n,' followed by
# as search for '(?<=\d)\s+(?=\d)' and replace with ','
# or be lazy and use this function
str_func <- function(x, n_digist = 16){
  tmp <- capture.output(print(c(x), digits = n_digist))
  tmp <- sapply(tmp, gsub, pattern = "\\s*\\[1\\]\\s", replacement = "", USE.NAMES = F)
  tmp <- sapply(tmp, gsub, pattern = "\\s*\\[\\d+\\]\\s", replacement = ",\\ ", USE.NAMES = F)
  tmp <- sapply(tmp, gsub, pattern = "(?<=\\d)\\s+(?=\\d|-)", replacement = ",\\ ", perl = T, USE.NAMES = F)

  tmp <- paste0(c("c(", tmp, " )"), collapse = "")

  max_lengt <- floor(8191 * .75)
  n_nums_before_break = floor(max_lengt / (n_digist + 4))
  gsub(pattern = paste0("((\\d,\\ .*?){", n_nums_before_break - 1, "})(,\\ )"), replacement = "\\1,\n\\ ",
       x = tmp, perl = T)
}

get_expect_equal <- function(x, eps){
  arg_name <- deparse(substitute(x))
  expects <- unlist(lapply(x, str_func))
  tol_string = if(!missing(eps)) paste0("\n, tolerance = " , eval(bquote(.(eps)))) else ""
  expects <- mapply(function(e, index_name)
    paste0("expect_equal(c(", arg_name, "$", index_name, "),\n", e,
           tol_string, ")", collapse = ""),
    e = expects, index_name = names(expects))

  out <- paste0(c("{", paste0(expects, collapse = "\n\n"), "}\n"), collapse = "\n")
  cat(out)
  invisible()
}
