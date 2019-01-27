# this is a very bad name...
approximator <- function(..., start, do_checks = FALSE){
  objs <- list(...)
  stopifnot(length(objs) > 1)

  is_mvn <- sapply(objs, "[[", "is_mvn")
  if(any(!is_mvn)){
    # find mode
    require(nloptr)
    eval_f <- function(x)
      -sum(sapply(objs, function(z) z$f(x)))
    eval_grad_f <- function(x)
      -rowSums(sapply(objs, function(z) z$deriv(x)))

    opt <- nloptr(start, eval_f = eval_f, eval_grad_f = eval_grad_f,
                  opts = list(algorithm = "NLOPT_LD_LBFGS", xtol_rel = 1e-8,
                              check_derivatives = do_checks, print_level = 1))

    if(opt$status < 0)
      stop("failed with code ", opt$status)

    val <- opt$solution

  } else
    val <- numeric(length(start))

  # compute concetration matrix
  neg_Ks <- lapply(objs, function(z) z$n_hessian(val))
  neg_K <- Reduce("+", neg_Ks)
  Sig <- solve(neg_K)

  k1 <- numeric(length(start))
  if(any(!is_mvn)){
    os <- objs[!is_mvn]
    objs <- objs[is_mvn]

    for(i in seq_along(os))
      k1 <- k1 + neg_Ks[[i]] %*% val + os[[i]]$deriv(val)

  }

  out <- function(...)
    {
      ds <- list(...)
      stopifnot(length(ds) == length(objs))
      for(i in seq_along(ds))
        k1 <- k1 + objs[[i]]$deriv_z(ds[[i]])

      mu = drop(solve(neg_K, k1))
      Sig = Sig

      require(mvtnorm)
      p_ldens <- dmvnorm
      formals(p_ldens)[c("mean", "sigma", "log")] <- list(mu, Sig, TRUE)
      p_sample <- rmvnorm
      formals(p_sample)[c("n", "mean", "sigma")] <- list(1L, mu, Sig)

      list(p_ldens = p_ldens, p_sample = p_sample, mu = mu, Sig = Sig)
    }
  structure(out, Sig = Sig, val = val)
}
