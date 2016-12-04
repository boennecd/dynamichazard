## ----setup, include=FALSE------------------------------------------------
knitr::knit_hooks$set(
  mySettings  = function(before, options, envir){
    if (before && options$mySettings){ 
      par(
        mar = c(10, 10, 4, 4),
        bty = "n",
        xaxs = "i",
        pch=16,
        cex= (cex <- .4),
        cex.axis = .8 / cex,
        cex.lab = .8 / cex,
        lwd= 1)
      options(digits = 3, width = 80, warn = -1)
    }})

knitr::opts_chunk$set(echo = TRUE, mySettings=TRUE, fig.height=3.5, fig.width = 6,
                      warning = F, message = F)
knitr::opts_knit$set(warning = F, message = F)

## ----echo=FALSE---------------------------------------------------------------
current_sha <- httr::content(
  httr::GET("https://api.github.com/repos/boennecd/dynamichazard/git/refs/heads/master")
  )$object$sha

stopifnot(length(current_sha) > 0 && class(current_sha) == "character")

current_version <- paste0("boennecd/dynamichazard@", current_sha)

## -----------------------------------------------------------------------------
current_version

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github(current_version)

