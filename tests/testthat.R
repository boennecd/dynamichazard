Sys.setenv("R_TESTS" = "")

library(testthat)
library(biglm)
library(dynamichazard)

options(ddhazard_max_threads = 1)
test_check("dynamichazard")
