binom <- function(y, X, offset, family){
  list(
    is_mvn = FALSE,
    f     = function(coefs)
    {
      mu <- family$linkinv(drop(X %*% coefs) + offset)
      sum(ifelse(y, log(mu), log(1 - mu)))
    },
    deriv = function(coefs)
    {
      eta <- drop(X %*% coefs + offset)
      mu <- family$linkinv(eta)
      d_f <- ifelse(
        y, 1 / (mu + .Machine$double.eps),
        - 1 / (1 - mu + .Machine$double.eps))
      colSums(X * d_f * family$mu.eta(eta))
    },
    n_hessian = function(coefs)
    {
      exp_eta <- exp(drop(X %*% coefs) + offset)
      exp_e_p1 <- 1 + exp_eta
      dd_f <- exp_eta / (exp_e_p1 * exp_e_p1 + .Machine$double.eps)
      crossprod(X * sqrt(dd_f))
    }
  )
}
