test_that("predict works for second order random walk",
          expect_true(FALSE))
test_that("residuals works for second order random walk",
          expect_true(FALSE))

test_that("Implement state space errors for EKF and exponential model. Requires computation of lag one correlation",
          expect_true(FALSE))

test_that("How to compute the starting value for lag one cov with UKF",
          expect_true(FALSE))

get_design_matrix <- function(...) environment(ddhazard)$get_design_matrix(...)

test_that("Implement lag-one-cov with weights", {
  expect_true(FALSE)
})

test_that("Works with exponential model for post_approx",
          expect_true(FALSE))

test_that("Gets same result regardless of which exponential model specification you give with post_approx",
          expect_true(FALSE))

test_that("This work when get_design_matrix when function is defined not in global scope", {
  # Simulate data
  set.seed(11111)
  sims = as.data.frame(test_sim_func_logit(n_series = 5e2, n_vars = 3, beta_start = 1,
                                           intercept_start = - 5, sds = c(sqrt(.1), rep(1, 3)),
                                           x_range = 1, x_mean = .5)$res)

  dum <- function(x) cbind(x, -x)
  design_with_fixed <- get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + ddFixed(dum(x2)) + x3), sims)

})


test_that("Figure out why Q gets big in this example", {
  set.seed(849239)
  sims_exp <- test_sim_func_exp(n_series = 4e2, n_vars = 3, x_range = 1, t_max = 10, x_mean = 0.5,
                                beta_start = 1, intercept_start = -4)

  fit <- suppressWarnings(ddhazard(form, Q_0 = diag(10, 2), Q = diag(1, 2),
                                   data = sims_exp$res, id = sims_exp$res$id,
                                   by = 1, model = "exp_bin", max_T = 10,
                                   control = list(fixed_terms_method = "E_step",
                                                  save_risk_set = F, save_data = F,
                                                  method = "UKF", debug = T), verbose = 5))
})
