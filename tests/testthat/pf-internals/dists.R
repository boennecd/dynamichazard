odist <- function(y, X, offset, family){
  if(inherits(family, "family") && family$family == "binomial" &&
     family$link == "logit"){
    linkinv <- family$linkinv
    l <- function(y, mu)
      sum(ifelse(y, log(mu), log(1 - mu)))

    mu.eta <- family$mu.eta
    lp <- function(y, eta){
      mu <- linkinv(eta)
      d_f <- ifelse(
        y, 1 / (mu + .Machine$double.eps),
        - 1 / (1 - mu + .Machine$double.eps))
      d_f * mu.eta(eta)
    }

    nllp <- function(y, eta){
      exp_eta <- exp(eta)
      exp_e_p1 <- 1 + exp_eta
      exp_eta / (exp_e_p1 * exp_e_p1 + .Machine$double.eps)
    }

  } else if(inherits(family, "family") && family$family == "binomial" &&
            family$link == "cloglog") {
    linkinv <- family$linkinv
    l <- function(y, mu)
      sum(ifelse(y, log(mu), log(1 - mu)))

    mu.eta <- family$mu.eta
    lp <- function(y, eta){
      mu <- linkinv(eta)
      d_f <- ifelse(
        y, 1 / (mu + .Machine$double.eps),
        - 1 / (1 - mu + .Machine$double.eps))
      d_f * mu.eta(eta)
    }

    nllp <- function(y, eta){
      exp_eta <- exp(eta)
      neg_exp_eta_expm1 <- -expm1(-exp_eta)
      -ifelse(
        y,
        exp(eta - exp_eta) * (1 - exp_eta) / neg_exp_eta_expm1 -
          exp(2 * eta - 2 * exp_eta) / neg_exp_eta_expm1^2 ,
        - exp_eta)
    }

  } else if(family == "exponential") {
    stopifnot(inherits(y, "Surv") && attr(y, "type") == "right")

    linkinv <- function(x)
      exp(-x)
    l <- function(y, mu)
      sum(ifelse(y[, 2], -log(mu) - y[, 1] / mu, - y[, 1] / mu))

    lp <- function(y, eta){
      ifelse(
        y[, 2],
        1 - exp(eta) * y[, 1],
        - exp(eta) * y[, 1])
    }

    nllp <- function(y, eta)
      exp(eta) * y[, 1]

  } else
    stop("not implemented")

  list(
    is_mvn = FALSE,
    f     = function(coefs)
    {
      mu <- linkinv(drop(X %*% coefs) + offset)
      l(y, mu)
    },
    deriv = function(coefs)
    {
      eta <- drop(X %*% coefs + offset)
      colSums(X * lp(y, eta))
    },
    n_hessian = function(coefs)
    {
      eta <- drop(X %*% coefs) + offset
      crossprod(X * nllp(y, eta), X)
    }
  )
}
