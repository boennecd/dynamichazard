Sys.setenv("R_TESTS" = "")

library(testthat)
library(biglm)
library(dynamichazard)

is_build_win = T
if(is_build_win)
  options(ddhazard_max_threads = 1)

test_check("dynamichazard")
