if(interactive()){
  library(testthat)
}

# Had issues with win builder. Thus, these lines
test_name <- "options"
cat("\nRunning", test_name, "\n")

test_that("'ddhazard_max_threads' defaults to -1",{
  old_opt <- options()
  options(ddhazard_max_threads = NULL)
  with(environment(ddhazard), .onLoad())
  expect_equal(getOption('ddhazard_max_threads'), -1)
  options(old_opt)
})



# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
