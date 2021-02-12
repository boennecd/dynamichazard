# Had the same issue as in this thread: https://github.com/hadley/testthat/issues/86
# or maybe https://github.com/r-lib/devtools/issues/1526
Sys.unsetenv("R_TESTS")

options(deparse.max.lines = 5,  testthat.summary.max_reports = 1000L)
testthat::test_check("dynamichazard", reporter = "summary")
