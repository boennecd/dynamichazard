context("Testing options")

test_that("'ddhazard_max_threads' defaults to -1",{
  old_opt <- options()
  options(ddhazard_max_threads = NULL)
  with(environment(ddhazard), .onLoad())
  expect_equal(getOption('ddhazard_max_threads'), -1)
  options(old_opt)
})
