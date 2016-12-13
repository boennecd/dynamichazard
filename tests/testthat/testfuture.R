test_that("predict works for second order random walk",
          expect_true(FALSE))
test_that("residuals works for second order random walk",
          expect_true(FALSE))

test_that("Implement state space errors for EKF and exponential model. Requires computation of lag one correlation",
          expect_true(FALSE))

test_that("How to compute the starting value for lag one cov with UKF",
          expect_true(FALSE))

get_design_matrix <- function(...) environment(ddhazard)$get_design_matrix(...)


test_that("This work when get_design_matrix when function is defined not in global scope", {
  # Simulate data
  set.seed(11111)
  sims = as.data.frame(test_sim_func_logit(n_series = 5e2, n_vars = 3, beta_start = 1,
                                           intercept_start = - 5, sds = c(sqrt(.1), rep(1, 3)),
                                           x_range = 1, x_mean = .5)$res)

  dum <- function(x) cbind(x, -x)
  design_with_fixed <- get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + ddFixed(dum(x2)) + x3), sims)

})
