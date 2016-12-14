test_that("Works with EKF and logit model",
          expect_true(F))

test_that("Works with UKF and logit model",
          expect_true(F))

test_that("Works with EKF and continous time model",
          expect_true(F))

test_that("Works with UKF and continous time model",
          expect_true(F))

test_that("Works with second order random walk and logit model",
          expect_true(F))

test_that("Works with second order random walk and continous time model",
          expect_true(F))

test_that("Works when only one is time varying",
          expect_true(F))

test_that("Works when one is fixed",
          expect_true(F))
