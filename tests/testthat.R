library(testthat)
library(biglm)
library(dynamichazard)

is_build_win = T
if(is_build_win)
  options(ddhazard_max_threads = 1)

cat("Running tests on:\n")
print(R.version)

tryCatch(
  test_check("dynamichazard"),
  error = function(e){
    print(e$call)
    print(e$message)
    print(.Traceback)
    stop(e)
  })
