context("Testing test_utils")

set.seed(101010)

test_that("Testing util functions to sim for test", {
  test_pairs <- cbind("new funcs" = c(get_exp_draw(), get_unif_draw(), get_norm_draw()),
                     "funcs" = c(rexp, runif, rnorm))

  for(i in seq_len(nrow(test_pairs))){
    new_func <- test_pairs[i, 1][[1]]
    func <- test_pairs[i, 2][[1]]

    expect_equal(length(new_func(1)), 1)
    expect_equal(length(new_func(10)), 10)
    expect_equal(length(new_func(100)), 100)

    seed <- 101
    set.seed(seed)
    new_func(re_draw = T)
    res1 <- new_func(10)
    set.seed(seed)
    res2 <- func(10)
    expect_equal(res1, res2)
  }
})

test_that("Testing util functions to sim series for tests", {
  skip_on_cran()

  set.seed(4321)
  n_series <- 1e4
  t_max <- 10

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


test_that("Whether covs are fixed or not depends on is fixed argument", {
  set.seed(1999293)
  inter_start <- rnorm(1)
  beta_start <- rnorm(10)
  for(func in c(test_sim_func_logit, test_sim_func_exp)){
    sims <- func(n_series = 1e1, n_vars = 10, beta_start = beta_start,
                 intercept_start = inter_start, t_max = 10)

    expect_true(all(apply(sims$betas, 2, function(x) sum(duplicated(x)) == 0)))

    sims <- test_sim_func_exp(n_series = 1e2, n_vars = 10, beta_start = beta_start,
                              intercept_start = inter_start, t_max = 10, is_fixed = c(1, 4))
    expect_equal(apply(sims$betas, 2, function(x) sum(duplicated(x)) == 0),
                 !seq_len(11) %in% c(1, 4))
  }
})

test_that("Sim functions gives previous results", {
  set.seed(897730)
  sim_logit <- test_sim_func_logit(1e2, n_vars = 3)
  sim_exp <- test_sim_func_exp(1e2, n_vars = 3)

  # save_to_test(sim_logit, "util_sim_logit")
  # save_to_test(sim_exp, "util_sim_exp")

  expect_equal(sim_logit, read_to_test("util_sim_logit"), tolerance = 1.490116e-08)
  expect_equal(sim_exp, read_to_test("util_sim_exp"), tolerance = 1.490116e-08)
})

test_that("`linear_mapper`s gives expected result", {
  set.seed(56219373)
  A <- structure(
    c(0.157, 0.025, 1.14, -2.055, -0.376, -0.837, -0.414,
      -0.627, 0.272, -0.574, 0.526, -0.088, 0.566, 0.809, -0.198),
    .Dim = c(5L, 3L))
  x <- rnorm(3)
  z <- rnorm(5)
  X <- matrix(rnorm(3 * 3), 3)
  Z <- matrix(rnorm(5 * 5), 5)

  test_expr <- quote({
    expect_equal(out$A, mat)
    expect_equal(out$A_x, mat %*% x)
    if(has_z)
      expect_equal(out$A_T_z, t(mat) %*% z)

    expect_equal(out$A_X, mat %*% X)
    expect_equal(out$X_A_T, X %*% t(mat))
    expect_equal(out$A_X_A_T, mat %*% X %*% t(mat))

    if(has_z) {
      expect_equal(out$A_T_Z, t(mat) %*% Z)
      expect_equal(out$Z_A, Z %*% mat)
      expect_equal(out$A_T_Z_A, t(mat) %*% Z %*% mat)
    }
  })

  out <- linear_mapper_test(
    A, x, X, z, Z, type = "dens_mapper", R = matrix()[0, 0, drop = FALSE])
  eval(do.call(substitute, list(
    test_expr, list(mat = quote(A), has_z = TRUE))))

  R <- diag(5)[, c(1, 4, 5)]
  out <- linear_mapper_test(
    R, x, X, z, Z, type = "select_mapper", R = matrix()[0, 0, drop = FALSE])
  eval(do.call(substitute, list(
    test_expr, list(mat = quote(R), has_z = TRUE))))

  C <- matrix(c(
     7,  2,  1,
     0,  3, -1,
    -3,  4, -2), byrow = TRUE, ncol = 3)
  out <- linear_mapper_test(
    C, x, X, x, t(X), type = "inv_mapper", R = matrix()[0, 0, drop = FALSE])
  eval(do.call(substitute, list(test_expr, list(
    z = quote(x), Z = t(X), mat = quote(solve(C)), has_z = TRUE))))

  R <- diag(3)[, c(1, 3)]
  out <- linear_mapper_test(
    C, x, X, x, t(X), type = "inv_sub_mapper", R = R)
  eval(do.call(substitute, list(test_expr, list(
    mat = quote(t(R) %*% solve(C)), has_z = FALSE))))
})
