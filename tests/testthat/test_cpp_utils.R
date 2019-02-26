context("Running test_cpp_utils.R:")

test_that("lambert_W0 gives close to same results as gsl::lambert_W0", {
  skip_if_not_installed("gsl")

  x <- log(seq(exp(1e-12), exp(3e-2), length.out = 1e3))

  expect_equal(sapply(x, lambert_W0_test), gsl::lambert_W0(x), tolerance = 5e-5)
})

test_that("trunc_lp_in_exponential_dist does not truncate when not needed", {
  skip_on_cran()
  .interactive <- interactive()

  eps <- exp(trunc_eta_exponential_test_log_eps())

  test_func <- function(eta, .t, is_event){
    q <- bquote({
      ans <- trunc_eta_exponential_test(
        eta = .(eta), at_risk_length = .(.t), is_event = .(is_event))

      if(.(.interactive)){
        cat(c(.(eta), .(.t), .(is_event)), "\n")
        cat("Likelihood is", exp(.(eta) * .(is_event) - exp(.(eta)) * .(.t)), "\n")
        cat("ans is ", unlist(ans), "\n")

      }
      should_trunc <- .(eta) * .(is_event) - exp(.(eta)) * .(.t) < log(eps)
      if(!should_trunc){
        expect_equal(ans$exp_eta_trunc, exp(.(eta)))
        expect_equal(ans$eta_trunc, .(eta))

      } else{
        expect_equal(exp(ans$eta_trunc), ans$exp_eta_trunc)
        expect_equal(ans$eta_trunc * .(is_event) - exp(ans$eta_trunc) * .(.t), log(eps))

        if(.(is_event)){
          expect_equal(
            ans$eta_trunc < - exp(ans$eta_trunc) * .(.t),
            .(eta) < - exp(.(eta)) * .(.t))
        }

      }
    })

    eval(q, envir = environment())
  }

  test_vals <- expand.grid(
    eta = -25:25,
    t = c(1e-2, 1, 1e2),
    is_event = c(TRUE, FALSE))

  invisible(
    mapply(test_func, eta = test_vals$eta, .t = test_vals$t, is_event = test_vals$is_event)
  )
})


test_that("round_if_almost_eq rounds to nearest boundary as expected", {
  set.seed(56219385)
  n <- 1e3

  x <- replicate(n, {
    x <- 1
    for(i in 1:10){
      tmp <- rexp(1, 1)
      n <- sample.int(1e2, 1)
      delta <- tmp / n
      for(j in 1:n)
      x <- x - delta
      x <- x + tmp
    }
    x
  })

  boundaries <- seq(0, 2, .1)
  x_ord <- order(x) - 1L
  x_trans <- drop(round_if_almost_eq(x, x_ord, boundaries = boundaries))

  correct_bound <- boundaries[11]
  expect_gt(sx <- sum(x_trans == correct_bound), sum(x == correct_bound))
  expect_equal(sx, length(x))

  # values not at the boundary should not get truncated
  x <- rnorm(n)
  x_ord <- order(x) - 1L
  x_trans <- drop(round_if_almost_eq(x, x_ord, boundaries = boundaries))
  expect_equal(x_trans, x)
})

test_that("rep_vec gives expected result", {
  n_i <- 1e2
  js <- rnorm(1e2)
  expect_equal(c(sapply(js, "*", rep.int(1, n_i))), drop(rep_vec(js, n_i)))

  # microbenchmark::microbenchmark(
  #   R0 <- c(sapply(js, "*", rep.int(1, n_i))),
  #   R1 <- c(sapply(js, rep, times = n_i)),
  #   R2 <- c(t.default(matrix(1L, length(js), n_i) * js)),
  #   R3 <- rep_vec(js, n_i))
  # all.equal(R0, R1)
  # all.equal(R1, R2)
  # all.equal(R2, R3)
})



test_that("selection_matrix_map gives correct results", {
  L <- matrix(0, 2, 3)
  L[2, 1] <- L[1, 3] <- 1

  x <- c(0.64, 0.41, 0.99)
  X <- structure(c(
    -0.08, -1.83, 1.99, -0.19, -1.55, 0.9, 0.1, 2.05, 0.04, 0.2, -0.03,
    -0.37, -0.76, -0.16, 1.44, -1.03, -0.36, 0.2, -0.03, 0.54, 0.15, -0.09,
    -0.18, 0.25), .Dim = c(3L, 8L))

  expect_equal(
    selection_matrix_map_vec_test(L, x, FALSE), L %*% x)
  expect_equal(
    selection_matrix_map_mat_test(L, X, FALSE, FALSE), L %*% X)
  expect_equal(
    selection_matrix_map_mat_test(L, t(X), TRUE, FALSE), t(X) %*% t(L))

  x <- x[1:2]
  X <- X[1:2, ]

  expect_equal(
    selection_matrix_map_vec_test(L, x, TRUE), t(L) %*% x)
  expect_equal(
    selection_matrix_map_mat_test(L, X, FALSE, TRUE), t(L) %*% X)
  expect_equal(
    selection_matrix_map_mat_test(L, t(X), TRUE, TRUE), t(X) %*% L)
})



test_that("'get_resample_idx_n_log_weight' gives correct results", {
  ws <- c(rep(1, 5), rep(2, 5))
  ws <- log(ws / sum(ws))

  rws <- 1:10
  rws <- log(rws / sum(rws))

  idx <- c(2, 2, 2, 4:6, 8, 8, 8, 8)

  count <- table(idx)
  idx_out <- as.integer(names(count)) - 1L
  ws_out <- ws[idx_out + 1L] - rws[idx_out + 1L] + log(count)
  ws_out <- exp(ws_out)
  ws_out <- log(ws_out / sum(ws_out))

  cpp_out <- test_get_resample_idx_n_log_weight(
    log_weights = ws, log_resample_weights = rws, resample_idx = idx - 1L)

  expect_equal(idx_out, drop(cpp_out$idx))
  expect_equal(c(ws_out), drop(cpp_out$log_weights), check.attributes = FALSE)
})
