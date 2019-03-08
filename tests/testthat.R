# Had the same issue as in this thread: https://github.com/hadley/testthat/issues/86
Sys.setenv("R_TESTS" = "")

options(deparse.max.lines = 5)
suppressWarnings(RNGversion("3.5.0"))

testthat::test_check("dynamichazard", reporter = "summary")
