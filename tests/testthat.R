Sys.setenv("R_TESTS" = "")

library(testthat)
library(biglm)
library(dynamichazard)

test_check("dynamichazard")
