parallelglm <- asNamespace("dynamichazard")$parallelglm
sim_expr <- expression({
  sim_func <- switch(
    family,
    binomial =function(eta)
      1 / (1 + exp(-eta)) >  runif(length(eta)),
    poisson = function(eta) rpois(n = n, exp(eta)))

  X <- structure(
    c(rep(1, n), runif(n*q, -1, 1)), dim = c(n, q + 1))
  .beta <- c(0, runif(q, -1, 1))
  eta <- drop(X %*% .beta)
  y <- sim_func(eta)

  offset <- if(use_offset)
    offset = runif(n, -1, 1) else rep(0, n)

  if(use_weights){
    weights <- runif(n)
    .weights <- weights / (sum(weights) / n)
  } else
    .weights <- rep(1, n)

  epsilon <- if(one_it) .Machine$double.xmax else 1e-8
})

n <- 1e5
q <- 25
family <- "poisson"
use_offset <- FALSE
use_weights <- FALSE
one_it <- FALSE
eval(sim_expr)

microbenchmark::microbenchmark(
  glm = coef_glm <- glm.fit(
    X, y, family = match.fun(family)(),
    weights = .weights, offset = offset,
    control = glm.control(epsilon = epsilon))$coefficients,

  parallelglm_chol = coef_parallelglm_Chol <- parallelglm(
    X = t(X), Ys = y,
    weights = .weights, offsets = offset, beta0 = numeric(),
    family = family, method = "quick",
    tol = epsilon, nthreads = 7),

  parallelglm_QR = coef_parallelglm_QR <- parallelglm(
    X = t(X), Ys = y,
    weights = .weights, offsets = offset, beta0 = numeric(),
    family = family, method = "QR",
    tol = epsilon, nthreads = 7),

  times = 3)

expect_equal(coef_glm, drop(coef_parallelglm_Chol))
expect_equal(coef_glm, drop(coef_parallelglm_QR))

#####
parallelglm_QR_get_R_n_f <-
  asNamespace("dynamichazard")$parallelglm_QR_get_R_n_f

n <- 1e4
q <- 2
family <- "poisson"
use_offset <- FALSE
use_weights <- FALSE
one_it <- FALSE
beta0 <- rep(0, q + 1)
eval(sim_expr)

args <- list(
  X = t(X), Ys = y, weights = .weights, offsets = offset, beta0 = beta0,
  family = family, nthreads = 7)


options(digits = 3)
microbenchmark::microbenchmark(
  do.call(parallelglm_QR_get_R_n_f, c(args, list(block_size =   500))),
  do.call(parallelglm_QR_get_R_n_f, c(args, list(block_size =  1000))),
  do.call(parallelglm_QR_get_R_n_f, c(args, list(block_size =  5000))),
  do.call(parallelglm_QR_get_R_n_f, c(args, list(block_size = 10000))),
  times = 25)
