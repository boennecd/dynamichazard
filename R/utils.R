.rbind_list <- function(l){
  out <- NULL
  for(i in 1:length(l[[1]]))
    out <- c(out, list(do.call(.rbind_list_inner, lapply(l, "[[", i))))

  out <- data.frame(out)
  names(out) <- names(l[[1]])
  out
}

.rbind_list_inner <- function(...){
  args <- list(...)
  if(is.factor(args[[1]]))
    # see https://stackoverflow.com/a/3449403/5861244
    return(factor(do.call(c, lapply(args, as.character))))

  do.call(c, args)
}
