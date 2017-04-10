cur_load = if(exists(".onLoad()")) .onLoad else function() { NULL }
.onLoad <- function(libname, pkgname){
  cur_load()

  op <- options()

  may_have_speed_glm <- FALSE
  try({
    find.package("speedglm")
    may_have_speed_glm <- TRUE
  })

  op.dynhazard <- list(
    ddhazard_max_threads = -1,
    ddhazard_use_speedglm = may_have_speed_glm)
  toset <- !(names(op.dynhazard) %in% names(op))
  if(any(toset)) options(op.dynhazard[toset])

  invisible()
}
