context("Testing cpp sampling functions")

test_that("sample_indices give the same as R with same seed", {
  probs <- 1:20 / sum(1:20)

  if(!exists(".Random.seed"))
    set.seed(NULL)
  for(i in 1:1e2){
    x <- .Random.seed
    R_res <- sample.int(length(probs), replace = TRUE, prob = probs)
    .new <- .Random.seed

    .Random.seed <<- x
    cpp_res <- sample_indices_test(length(probs), probs)

    expect_equal(R_res - 1, drop(cpp_res))
    expect_equal(.new, .Random.seed)
  }

  for(i in 1:1e2){
    x <- .Random.seed
    R_res <- sample.int(length(probs), size = 2 * length(probs),replace = TRUE, prob = probs)
    .new <- .Random.seed

    .Random.seed <<- x
    cpp_res <- sample_indices_test(2 * length(probs), probs)

    expect_equal(R_res - 1, drop(cpp_res))
    expect_equal(.new, .Random.seed)
  }
})

test_that("mvrnorm gives the same results with same seed", {
  n <- 5
  mu <- -2:2
  Sigma <- matrix(.5, ncol=n, nrow = n)
  diag(Sigma) <- 1:n

  x <- .Random.seed
  r1 <- mvrnorm_test(mu = mu, sigma_chol = chol(Sigma))
  .new <- .Random.seed
  .Random.seed <<- x
  r2 <- mvrnorm_test(mu = mu, sigma_chol = chol(Sigma))

  expect_equal(r1, r2)
  expect_equal(.new, .Random.seed)
})

test_that("mvrnorm gives expected sample mean and variance", {
  skip_on_cran()
  n <- 3
  mu <- -1:1

  set.seed(83685361)
  for(i in 1:10){
    A <- matrix(runif(n^2)*2-1, ncol=n)
    Sigma <- t(A) %*% A

    samp <- t(replicate(1000, drop(
      mvrnorm_test(mu = mu, sigma_chol = chol(Sigma)))))

    expect_equal(colMeans(samp), mu, tolerance = .1)
    expect_equal(cov(samp), Sigma, tolerance = .25)
  }
})

test_that("different seed gives different results", {
  n <- 5
  mu <- -2:2
  Sigma <- matrix(.5, ncol=n, nrow = n)
  diag(Sigma) <- 1:n

  set.seed(1)
  r1 <- mvrnorm_test(mu = mu, sigma_chol = chol(Sigma))
  set.seed(2)
  r2 <- mvrnorm_test(mu = mu, sigma_chol = chol(Sigma))

  expect_false(isTRUE(all.equal(r1, r2)))
})

test_that("cpp systematic_resampling function gives the same as R version", {
  test_func <- function(n, probs){
    n_probs <- length(probs)
    u <- runif(1, 0, 1/n)
    u <- c(u, u + ((2:n) - 1) / n)
    p_cum_next <- cumsum(probs)
    p_cum <- c(0, head(p_cum_next, -1))
    sapply(
      u, function(x) which(p_cum <= x & x < p_cum_next))
  }

  for(n in c(10, 100, 1000)) {
    for(dum in 1:10){

      probs <- runif(100, 0, 1); probs <- probs / sum(probs)
      seed <- .Random.seed
      out <- test_func(n, probs)
      seed_after <- .Random.seed

      .Random.seed <<- seed
      out_cpp <- systematic_resampling_test(n, probs)

      expect_equal(seed_after, .Random.seed)
      expect_equal(out - 1, drop(out_cpp))

    }
  }
})

test_that("'mvtrnorm' yields close to uniform sample", {
  skip_on_cran()
  skip_if_not_installed("mvtnorm")

  n <- 2
  mu <- c(-1, 1)
  Sigma <- matrix(.5, ncol=n, nrow = n)
  diag(Sigma) <- rep(1, n)

  nu <- 5
  set.seed(6922317)
  X <- drop(replicate(
    1000, mvtrnorm_test(mu = mu, sigma_chol = chol(Sigma), nu = nu)))
  # plot(X[1, ], X[2, ], xlim = c(-6, 6), ylim = c(-6, 6))

  library(mvtnorm)
  p1 <- apply(X, 2, function(x)
    pmvt(lower = -Inf, upper = c(x[1], Inf), df = nu, corr = Sigma,
         type = "shifted", delta = mu))
  expect_true(ks.test(p1, "punif", 0, 1)$p.value > .1)

  p2 <- apply(X, 2, function(x)
    pmvt(lower = -Inf, upper = c(Inf, x[2]), df = nu, corr = Sigma,
         type = "shifted", delta = mu))
  expect_true(ks.test(p2, "punif", 0, 1)$p.value > .1)

  # try with e.g., the following instead
  #   X <- t(matrix(rnorm(1000 * 2), ncol = 2) %*% chol(Sigma)) + mu
})

test_that("cpp sample_n_count_replicas gives the right result", {
  ps <- exp(1:1000 / 1000)
  ps <- ps / sum(ps)

  set.seed(seed <- 101)
  x <- sample_indices_test(size = 100, probs = ps)
  x <- table(drop(x))

  set.seed(seed)
  y <- sample_n_count_replicas_indices_test(size = 100, probs = ps)

  expect_equal(as.numeric(names(x)),  y[, 1L])
  expect_equal(unclass(x),  y[, 2L], check.attributes = FALSE)

  set.seed(seed)
  x <- systematic_resampling_test(size = 100, probs = ps)
  x <- table(drop(x))

  set.seed(seed)
  y <- sample_n_count_replicas_systematic_test(size = 100, probs = ps)

  expect_equal(as.numeric(names(x)),  y[, 1L])
  expect_equal(unclass(x),  y[, 2L], check.attributes = FALSE)
})
