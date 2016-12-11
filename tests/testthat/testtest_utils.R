if(interactive()){
  rm(list = ls())
  source("./R/test_utils.R")
  library(testthat)
}


# Had issues with win builder. Thus, these lines
test_name <- "test_utils"
cat("\nRunning", test_name, "\n")


set.seed(101010)

test_that("Testing util functions to sim for test", {
  test_pairs <- cbind("new funcs" = c(get_exp_draw, get_unif_draw, get_norm_draw),
                     "funcs" = c(rexp, runif, rnorm))

  for(i in seq_len(nrow(test_pairs))){
    new_func <- test_pairs[i, 1][[1]]
    func <- test_pairs[i, 2][[1]]

    tmp <- unlist(replicate(10^3, new_func(1000)))
    expect_equal(length(unique(tmp)), length(tmp))

    expect_equal(length(new_func(1)), 1)
    expect_equal(length(new_func(10)), 10)

    seed <- 101
    set.seed(seed)
    new_func(re_draw = T)
    res1 <- new_func(10)
    set.seed(seed)
    res2 <- func(10)
    expect_equal(res1, res2)
  }
})

set.seed(4321)
n_series <- 1e4
t_max <- 10

test_that("Testing util functions to sim series for tests", {
  for(func in c(test_sim_func_logit, test_sim_func_exp)){
    tmp = func(n_series, t_max = t_max)$res

    expect_equal(unique(tmp[, "id"]), seq_len(n_series))

    expect_true(all(tapply(tmp[, "tstop"], tmp[, "id"], function(ts) sum(ts > 10)) < 2))
    expect_true(any(tapply(tmp[, "tstop"], tmp[, "id"], function(ts) sum(ts >= 10)) > 0))

    expect_true(all(tapply(tmp[, "event"], tmp[, "id"], function(ev) sum(ev)) < 2))
    expect_true(any(tapply(tmp[, "event"], tmp[, "id"], function(ev) sum(ev)) > 0))

    expect_true(all(tmp[, "tstart"] < tmp[, "tstop"]))

    expect_true(all(tmp[-1, "id"] != tmp[-nrow(tmp), "id"] |
                      tmp[-1, "tstart"] == tmp[-nrow(tmp), "tstop"]))

    expect_true(all(tmp[-1, "id"] == tmp[-nrow(tmp), "id"] |
                      tmp[-nrow(tmp), "tstop"] > 4 |
                      tmp[-nrow(tmp), "event"]))

    expect_true(all(!(tmp$res$event == 1 &
                        tmp$res$tstop > t_max)))
  }
})

test_that("The head_neck_cancer must be defined for the ddhazard tests",
          expect_true(exists("head_neck_cancer")))

test_that("The pbc2 must be defined for the ddhazard tests",
          expect_true(exists("pbc2")))


test_that("Covs are fixed or not depends on is fixed argument", {
  set.seed(1999293)
  inter_start <- rnorm(1)
  beta_start <- rnorm(10)
  for(func in c(test_sim_func_logit, test_sim_func_exp)){
    sims <- func(n_series = 1e2, n_vars = 10, beta_start = beta_start,
                 intercept_start = inter_start, t_max = 10)

    expect_true(all(apply(sims$betas, 2, function(x) sum(duplicated(x)) == 0)))

    sims <- test_sim_func_exp(n_series = 1e2, n_vars = 10, beta_start = beta_start,
                              intercept_start = inter_start, t_max = 10, is_fixed = c(1, 4))
    expect_equal(apply(sims$betas, 2, function(x) sum(duplicated(x)) == 0),
                 !seq_len(11) %in% c(1, 4))
  }
})




# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
