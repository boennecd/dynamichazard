#' @importFrom Rcpp loadModule Module
loadModule("dd_exponential")
loadModule("dd_logistic")

cur_load = if(exists(".onLoad()")) .onLoad else function() { NULL }
.onLoad <- function(libname, pkgname){
  cur_load()

  op <- options()

  # TODO: avoid this by doing as e.g., in `boot::boot` by using `getOption`
  op.dynhazard <- list(
    ddhazard_max_threads = -1)
  toset <- !(names(op.dynhazard) %in% names(op))
  if(any(toset))
    options(op.dynhazard[toset])

  invisible()
}
