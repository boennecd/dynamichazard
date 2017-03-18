# Had the same issue as in this thread: https://github.com/hadley/testthat/issues/86
Sys.setenv("R_TESTS" = "")

library(testthat)
library(biglm)
library(dynamichazard)

is_build_win = T
if(is_build_win)
  options(ddhazard_max_threads = 2)

cat("Running tests on:\n")
print(R.version)


test_check("dynamichazard")
