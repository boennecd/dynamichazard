context("Testing parallelglm vs. glm")

# irls_qrnewton =
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

n <- 1e5
q <- 3
X <- structure(
  c(rep(1, n), rnorm(n*q)), dim = c(n, q + 1))
beta <- c(0, runif(q, -1, 1))
eta <- drop(X %*% beta)

y <- 1 / (1 + exp(-eta)) >  runif(n)
glm_out <- glm.fit(X, y, family = binomial(), control = glm.control(trace = TRUE))
glm_out$coefficients
beta

parallelglm <- asNamespace("dynamichazard")$parallelglm
sink("tmp.txt")
out <- parallelglm(X = t(X), Ys = y,
                   weights = numeric(), offsets = numeric(), beta0 = numeric(),
                   family = "binomial",
                   tol = 1e-8)
sink()

mean((out - beta)^2)
mean((glm_out$coefficients - beta)^2)
rbind(
  drop(out),
  glm_out$coefficients)

# speedglm::speedglm.wfit(y, X, family = binomial())$coefficients
#
# sink("tmptmp.txt")
# drop(irls_qrnewton(X, y)$coefficients)
# sink()

microbenchmark::microbenchmark(
  glm.fit(X, y, family = binomial()),
  parallelglm(X = t(X), beta0 = rep(0, q + 1), Ys = y, weights = rep(1, n), family = "binomial",
              offsets = rep(0, n), nthreads = 7),
  times = 10)
