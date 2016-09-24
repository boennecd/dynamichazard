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
n_series <- 1e5
t_max <- 10

test_that("Testing util functions to sim series for tests", {
  for(func in c(test_sim_func_logit, test_sim_func_exp)){
    tmp = func(n_series, t_max = t_max)$res

    expect_equal(unique(tmp[, "id"]), seq_len(n_series))

    expect_true(all(tapply(tmp[, "tstop"], tmp[, "id"], function(ts) sum(ts > 10)) < 2))
    expect_true(any(tapply(tmp[, "tstop"], tmp[, "id"], function(ts) sum(ts > 10)) > 0))

    expect_true(all(tapply(tmp[, "event"], tmp[, "id"], function(ev) sum(ev)) < 2))
    expect_true(any(tapply(tmp[, "event"], tmp[, "id"], function(ev) sum(ev)) > 0))

    expect_true(all(tmp[, "tstart"] < tmp[, "tstop"]))

    expect_true(all(tmp[-1, "id"] != tmp[-nrow(tmp), "id"] |
                      tmp[-1, "tstart"] == tmp[-nrow(tmp), "tstop"]))

    expect_true(all(tmp[-1, "id"] == tmp[-nrow(tmp), "id"] |
                      tmp[-nrow(tmp), "tstop"] > 4 |
                      tmp[-nrow(tmp), "event"]))
  }
})
