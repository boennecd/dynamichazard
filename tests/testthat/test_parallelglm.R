context("Testing parallelglm vs. glm")

test_that("glm and parallelglm gives the same", {
  set.seed(4611691)
  parallelglm <- asNamespace("dynamichazard")$parallelglm

  grid_vals <- expand.grid(
    family = c("binomial", "poisson"),
    use_offset = c(TRUE, FALSE),
    use_weights = c(TRUE, FALSE),
    one_it = c(TRUE, FALSE))
  grid_vals$family <- as.character(grid_vals$family)

  test_expr <- expression({
    for(i in 1:nrow(grid_vals)){
      vals <- grid_vals[i, ]
      with(vals, {
          sim_func <- switch(
            family,
            binomial =function(eta)
              1 / (1 + exp(-eta)) >  runif(length(eta)),
            poisson = function(eta) rpois(n = n, exp(eta)))

          X <<- structure(
            c(rep(1, n), runif(n*q, -1, 1)), dim = c(n, q + 1))
          .beta <<- c(0, runif(q, -1, 1))
          eta <<- drop(X %*% .beta)
          y <<- sim_func(eta)

          offset <- if(use_offset)
            offset = runif(n, -1, 1) else rep(0, n)

          if(use_weights){
            weights <- runif(n)
            .weights <<- weights / (sum(weights) / n)
          } else
            .weights <<- rep(1, n)

          epsilon <- if(one_it) .Machine$double.xmax else 1e-8

          glm_out <<- suppressWarnings( # glm gives warnings with non-integer weights
            glm.fit(
              X, y, family = match.fun(family)(),
              weights = .weights, offset = offset,
              control = glm.control(epsilon = epsilon))$coefficients)

          out <<- parallelglm(X = t(X), Ys = y,
                              weights = .weights, offsets = offset, beta0 = numeric(),
                              family = family,
                              tol = epsilon, nthreads = getOption("ddhazard_max_threads"))
      })
      expect_equal(glm_out, drop(out), tolerance = 1e-5)
    }})

  n <- 200; q <- 3
  eval(test_expr)

  n <- 2e3
  eval(test_expr)

  skip_on_cran()

  n <- 1e5; q <- 10
  eval(test_expr)
})


# irls_newton =
#   function(A, b, family=binomial, maxit=25, tol=1e-08)
#   {
#     s  = rep(0, ncol(A))
#     for(j in 1:maxit)
#     {
#       t      = A %*% s
#       g      = family()$linkinv(t)
#       gprime = family()$mu.eta(t)
#       z      = t + (b - g) / gprime
#       W      = as.vector(gprime^2 / family()$variance(g))
#
#       print(paste("eta = ", t,  "g = ", g, " gprime = ", gprime, " z = ", z, " W = ",  W))
#
#       wmin   = min(W)
#       if(wmin < sqrt(.Machine$double.eps))
#         warning("Tiny weights encountered")
#       s_old   = s
#
#       print(base::t(A) %*% (W * A))
#       print(crossprod(A,W*z))
#
#       s   = solve(base::t(A) %*% (W * A), crossprod(A,W*z))
#
#       print(s)
#
#       if(sqrt(crossprod(s - s_old)) < tol) break
#     }
#     list(coefficients=s,iterations=j)
#   }
#
# irls_qrnewton
