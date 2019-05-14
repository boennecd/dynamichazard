# Had the same issue as in this thread: https://github.com/hadley/testthat/issues/86
Sys.setenv("R_TESTS" = "")

options(deparse.max.lines = 5,  testthat.summary.max_reports = 1000L)
testthat::test_check("dynamichazard", reporter = "summary")
