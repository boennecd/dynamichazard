context("Testing parallelglm vs. glm")

# assign simulation function for later
get_sims <- function(family, use_offset, use_weights, n, q){
  sim_func <- switch(
    family,
    binomial = function(eta)
      1 / (1 + exp(-eta)) >  runif(length(eta)),
    poisson = function(eta) rpois(n = n, exp(eta)),
    cloglog = function(eta) - expm1(-exp(eta)) > runif(length(eta)))

  X <- structure(
    c(rep(1, n), runif(n * q, -1, 1)), dim = c(n, q + 1))
  .beta <- c(0, runif(q, -1, 1))
  eta <- drop(X %*% .beta)
  y <- sim_func(eta)

  offset <- if(use_offset)
    offset = runif(n, -1, 1) else rep(0, n)

  weights. <- if(use_weights){
    weights. <- runif(n)
    weights. / (sum(weights.) / n)
  } else
    rep(1, n)

  list(X = X, y = y, offset = offset, weights. = weights., eta = eta)
}

test_that("glm and parallelglm gives the same", {
  set.seed(4611691)
  grid_vals <- expand.grid(
    method = c("quick", "QR"),
    family = c("binomial", "poisson", "cloglog"),
    use_offset = c(TRUE, FALSE),
    use_weights = c(TRUE, FALSE),
    one_it = c(TRUE, FALSE),
    stringsAsFactors = FALSE)

  test_expr <- expression({
    for(i in 1:nrow(grid_vals)){
      vals <- grid_vals[i, ]
      sim <- with(
        vals, get_sims(family, use_offset, use_weights, n, q))

      epsilon <- if(vals$one_it) .Machine$double.xmax else 1e-8

      fam <- switch(
        vals$family, binomial = binomial(), cloglog = binomial("cloglog"),
        poisson = poisson())

      glm_out <- suppressWarnings( # glm gives warnings with non-integer weights
        glm.fit(
          sim$X, sim$y, family = fam,
          weights = sim$weights., offset = sim$offset,
          control = glm.control(epsilon = epsilon)))

      if(!vals$one_it)
        # force one more iteration. I only evaluate the deviance after a fit
        glm_out <- suppressWarnings(
          glm.fit(
            sim$X, sim$y, family = fam, start = glm_out$coefficients,
            weights = sim$weights., offset = sim$offset,
            control = glm.control(epsilon = epsilon)))

      out <- parallelglm(
        X = t(sim$X), Ys = sim$y, weights = sim$weights., offsets = sim$offset,
        beta0 = numeric(), family = vals$family, tol = epsilon,
        nthreads = getOption("ddhazard_max_threads"), method = vals$method,
        it_max = 100)

      expect_equal(glm_out$coefficients, c(out), tolerance =
                     if(vals$method == "quick") 1e-4 else 1e-5)
    }})

  n <- 200; q <- 3
  eval(test_expr)
  vals

  skip_on_cran()

  n <- 2e3
  eval(test_expr)

  n <- 1e5; q <- 10
  eval(test_expr)
})
