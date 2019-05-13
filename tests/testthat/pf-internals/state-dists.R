get_Q_0 <- function(Qmat, Fmat){
  eg  <- eigen(Fmat)
  las <- eg$values
  if(any(Mod(las) >= 1))
    stop("Divergent series")
  U   <- eg$vectors
  T. <- solve(U, t(solve(U, Qmat)))
  Z   <- T. / (1 - tcrossprod(las))
  out <- tcrossprod(U %*% Z, U)
  if(is.complex(out)){
    if(all(abs(Im(out)) < .Machine$double.eps^(3/4)))
      return(Re(out))

    stop("Q_0 has imaginary part")
  }

  out
}

fw <- function(parent, F., Q){
  n_hes <- solve(Q)
  QiF <- solve(Q, F.)
  list(
    is_mvn = TRUE,
    f         = function(coefs)
    {
      require(mvtnorm)
      dmvnorm(coefs, F. %*% parent, Q, log = TRUE)
    },
    deriv     = function(coefs) -solve(Q, coefs - F. %*% parent),
    deriv_z   = function(p) QiF %*% p,
    n_hessian = function(...) n_hes)
}

bw <- function(child, F., Q){
  n_hes <- crossprod(F., solve(Q, F.))
  FtQi <- crossprod(F., solve(Q))
  list(
    is_mvn = TRUE,
    f         = function(coefs)
    {
      require(mvtnorm)
      dmvnorm(child, F. %*% coefs, Q, log = TRUE)
    },
    deriv     = function(coefs) FtQi %*% (child - F. %*% coefs),
    deriv_z   = function(ci) FtQi %*% ci,
    n_hessian = function(...) n_hes)
}

prior <- function(F., Q, Q_0, mu_0){
  if(missing(Q_0))
    Q_0 <- get_Q_0(Qmat = Qmat, Fmat = Fmat)
  if(missing(mu_0))
    mu_0 <- numeric(nrow(F.))
  Qs_mus  <- list(list(Q = Q_0, mu = mu_0))

  get_mu_n_Q <- function(ti){
    stopifnot(is.integer(ti), ti > 0L)
    tip <- as.integer(ti) + 1L
    n_ele <- length(Qs_mus)
    if(length(Qs_mus) < tip){
      for(i in (n_ele + 1L):tip){
        old <- Qs_mus[[i - 1L]]
        Q_new <- tcrossprod(F. %*% old$Q, F.) + Q
        Qs_mus[[i]] <<- list(Q = Q_new, mu = F. %*% old$mu)
      }
    }

    Qs_mus[[tip]]
  }

  function(ti){
    tmp <- get_mu_n_Q(ti)
    Qt <- tmp$Q
    mut <- tmp$mu

    n_hes <- solve(Qt)
    dz <- solve(Qt, mut)
    list(
      is_mvn = TRUE,
      f         = function(coefs)
      {
        require(mvtnorm)
        dmvnorm(coefs, mut, Qt, log = TRUE)
      },
      deriv     = function(coefs) -solve(Qt, coefs - mut),
      deriv_z   = function(...) dz,
      n_hessian = function(...) n_hes)
  }
}
