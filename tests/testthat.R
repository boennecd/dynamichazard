# Had the same issue as in this thread: https://github.com/hadley/testthat/issues/86
Sys.setenv("R_TESTS" = "")

testthat::test_check("dynamichazard")
