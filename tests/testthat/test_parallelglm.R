context("Testing parallelglm vs. glm")

# assign simulation function for later
get_sims <- function(family, use_offset, use_weights, n, q){
  sim_func <- switch(
    family,
    binomial = function(eta)
      1 / (1 + exp(-eta)) >  runif(length(eta)),
    poisson = function(eta) rpois(n = n, exp(eta)))

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

  list(X = X, y = y, offset = offset, weights. = weights.)
}

test_that("parallel QR update gives same cofficiens results with one iteration", {
  grid_vals <- expand.grid(
    family = c("binomial", "poisson"),
    use_offset = c(TRUE, FALSE),
    use_weights = c(TRUE, FALSE),
    stringsAsFactors = FALSE)

  beta0 <- rep(0, 3)

  for(i in 1:nrow(grid_vals)){
    vals <- grid_vals[i, ]
    sims <- with(
      vals, get_sims(family, use_offset, use_weights, 200, 2))

    tmp <- with(sims, parallelglm_QR_get_R_n_f(
      t(X), y, vals$family, beta0, weights., offset,
      tol = .Machine$double.xmax, # one iteration,
      block_size = 100, nthreads = 1))
    piv <- drop(tmp$pivot)
    piv[piv + 1] <- seq_along(piv)
    beta_par <- lm.fit(tmp$R[, piv, drop = FALSE], tmp$f)$coefficients

    glm_out <- with(
      sims, {
        fa <- get(vals$family)

        etastart <- drop(X %*% beta0) + offset
        mustart <- fa()$linkinv(etastart)

        suppressWarnings(
          glm.fit(X, y, family = fa(), etastart = etastart,
                  mustart = mustart,
                  control = glm.control(maxit = 1), offset = offset,
                  weights = weights.))
      })

    expect_equal(glm_out$coefficients, beta_par, check.attributes = FALSE,
                 label = paste(c(vals), collapse = " "))
  }
})

test_that("glm and parallelglm gives the same", {
  set.seed(4611691)
  grid_vals <- expand.grid(
    method = c("quick", "QR"),
    family = c("binomial", "poisson"),
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

      glm_out <- suppressWarnings( # glm gives warnings with non-integer weights
        glm.fit(
          sim$X, sim$y, family = match.fun(vals$family)(),
          weights = sim$weights., offset = sim$offset,
          control = glm.control(epsilon = epsilon))$coefficients)

      out <- parallelglm(
        X = t(sim$X), Ys = sim$y, weights = sim$weights., offsets = sim$offset,
        beta0 = numeric(), family = vals$family, tol = epsilon,
        nthreads = getOption("ddhazard_max_threads"), method = vals$method)

      expect_equal(glm_out, drop(out), tolerance = 1e-5)
    }})

  n <- 200; q <- 3
  eval(test_expr)

  skip_on_cran()

  n <- 2e3
  eval(test_expr)

  n <- 1e5; q <- 10
  eval(test_expr)
})



