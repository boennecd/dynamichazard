if(interactive()){
  library(biglm)
  bigglm_updateQR_rcpp <- function(...)
    with(environment(ddhazard), bigglm_updateQR_rcpp(...))
  bigglm_regcf_rcpp <- function(...)
    with(environment(ddhazard), bigglm_regcf_rcpp(...))
  bigqr.init <- with(environment(bigglm), bigqr.init)
}

biglm_func <- function(formula, data, model = "logit", maxit=8, tolerance=1e-12){
  tt<-terms(formula)
  # beta <- start

  for (i in 1:maxit){
    firstchunk <- TRUE
    deviance<-0
    rss<-0
    data(reset=TRUE) # Returns the next chunk of data
    n<-0
    while(!is.null(chunk<-data(reset=FALSE))){
      n<-n+nrow(chunk) # how many rows we are using

      ##########
      ## Find X

      mf<-model.frame(tt,chunk)
      mm<-model.matrix(tt,mf) # Get the terms
      p<-NCOL(mm)
      # w <- rep(1,nrow(mm))

      ##########
      ## Init dummy entries if needed

      if (firstchunk) {
        qr<-bigqr.init(p) # simple function that create a list with dummy entries
        attach(qr)

        if(i == 1)
          beta <- rep(0, p)
      }
      #if (!identical(assn, attr(mm,"assign")))
      #    stop("model matrices incompatible")

      ##########
      ## Find y and offsets

      y<-model.response(mf)
      if(is.null(off<-model.offset(mf))) off<-0

      ##########
      ## Find linear predcitor etc.

      eta<- beta %*% mm

      ##########
      ## Update QR

      bigglm_updateQR_rcpp(D, rbar, thetab, ss, checked, tol,
                           model,
                           mm, eta, offset, y)

      ##update.bigqr:
      ##function (bigQR, X, y, w = NULL, singcheck = FALSE, add.intercept = FALSE)
      ##{
      ##    if (NCOL(X) + add.intercept != length(bigQR$D))
      ##        stop("Wrong number of columns")
      ##    if (length(y) != NROW(X))
      ##        stop("Wrong number of rows")
      ##    if (length(w) == 0)
      ##        w <- rep(1, length(y))
      ##    if (length(y) != length(w))
      ##        stop("`weights' has wrong length")
      ##    storage.mode(X) <- "double"
      ##    storage.mode(y) <- "double"
      ##    storage.mode(w) <- "double"
      ##    bigQR <- .Call("updateQR", X, y, w, bigQR, add.intercept)
      ##    if (singcheck)
      ##        bigQR <- .Call("singcheckQR", bigQR)
      ##    bigQR
      ##}

      ## So this boils down to call to the C code:
      ## updateQR(X, y, w, bigQR, add.intercept)

      #if(!is.null(beta)){ # <-- not needed
      #    deviance<-deviance+sum(family$dev.resids(y,mu,w))
      #    rss<-rss+sum((y-mu)^2/(w*family$variance(mu)))*(sum(w)/length(w))
      #if (sandwich){
      #    xx<-matrix(nrow=nrow(mm), ncol=p*(p+1))
      #    xx[,1:p]<-mm*(drop(z)-off)
      #    for(i in 1:p)
      #        xx[,p*i+(1:p)]<-mm*mm[,i]
      #    xyqr<-update(xyqr,xx,rep(0,nrow(mm)),ww*ww)
      #}
      #}
      firstchunk <- FALSE
    }
    # iwlm <- list(call=sys.call(-1), qr=qr, iterations=i, # Create "biglm" class
    #              assign=attr(mm,"assign"), terms=tt, converged=FALSE,
    #              n=n,names=colnames(mm), weights=weights,rss=rss)
    # if(sandwich)
    #   iwlm$sandwich <- list(xy=xyqr)
    # class(iwlm) <- "biglm"

    betaold <- beta
    beta <- bigglm_regcf_rcpp(D, rbar, thetab,
                              ss, checked, tol)

    ##coef.biglm
    ##function (object, ...)
    ##{
    ##    if (!object$qr$checked)
    ##        object$qr <- singcheck.bigqr(object$qr)
    ##    rval <- coef(object$qr)
    ##    rval[object$qr$D == 0] <- NA
    ##    names(rval) <- object$names
    ##    rval
    ##}

    ##coef.bigqr
    ##function (bigQR, nvar = NULL, ...)
    ##{
    ##    p <- length(bigQR$D)
    ##    if (is.null(nvar))
    ##        nvar <- p
    ##    if (nvar < 1 | nvar > p)
    ##        stop("Invalid value of `nvar'")
    ##    if (!bigQR$checked)
    ##        bigQR <- singcheck.bigqr(bigQR)
    ##    tmp <- .Fortran("regcf", as.integer(p), as.integer(p * p/2),
    ##        bigQR$D, bigQR$rbar, bigQR$thetab, bigQR$tol, beta = numeric(p),
    ##        nreq = as.integer(nvar), ier = integer(1), DUP = FALSE)
    ##    if (tmp$ier != 0)
    ##        stop("Error in REGCF: can't happen")
    ##    tmp$beta
    ##}

    ## So this boils down to a call to the Fortran code
    ## regcf(...)

    if (i >= maxit){
      if (!quiet) warning("ran out of iterations and failed to converge")
      break
    }

    if (!is.null(betaold)){
      delta <- betaold-beta
      if (max(abs(delta)) < tolerance){
        iwlm$converged<-TRUE
        break
      }
    }
  }

  return(beta)
}
