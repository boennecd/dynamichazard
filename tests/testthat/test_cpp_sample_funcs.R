context("Testing cpp sampling functions")

sample_indices <- asNamespace("dynamichazard")$sample_indices
mvrnorm <- asNamespace("dynamichazard")$mvrnorm

test_that("sample_indices give the same as R with same seed", {
  probs <- 1:20 / sum(1:20)

  if(!exists(".Random.seed"))
    set.seed(NULL)
  for(i in 1:1e2){
    x <- .Random.seed
    R_res <- sample.int(length(probs), replace = TRUE, prob = probs)

    .Random.seed <<- x
    cpp_res <- sample_indices(probs)

    expect_equal(R_res - 1, drop(cpp_res))
  }
})

test_that("mvrnorm gives the same results with same seed", {
  n <- 5
  mu <- -2:2
  A <- matrix(runif(n^2)*2-1, ncol=n)
  Sigma <- t(A) %*% A

  x <- .Random.seed
  r1 <- mvrnorm(mu = mu, sigma_chol = chol(Sigma))
  .Random.seed <<- x
  r2 <- mvrnorm(mu = mu, sigma_chol = chol(Sigma))

  expect_equal(r1, r2)
})

test_that("mvrnorm gives expected sample mean and variance", {
  n <- 3
  mu <- -1:1

  set.seed(83685361)
  for(i in 1:10){
    A <- matrix(runif(n^2)*2-1, ncol=n)
    Sigma <- t(A) %*% A

    samp <- t(replicate(1000, drop(mvrnorm(mu = mu, sigma_chol = chol(Sigma)))))

    expect_equal(colMeans(samp), mu, tolerance = .1)
    expect_equal(cov(samp), Sigma, tolerance = .25)
  }
})
