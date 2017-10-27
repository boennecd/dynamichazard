context("Testing shiny app")

test_that("Shiny app starts when no arguments are passed", {
  skip_if_not_installed("shiny")

  ddhazard_app(quietly = T)
})

test_that("Shiny throws error when invalid arguments are passed", {
  skip_if_not_installed("shiny")

  expect_error(ddhazard_app(invalid_arg = T),
               "These input arguments are not recognized: invalid_arg")
})

test_that("Shiny app starts when default arguments are passed", {
  skip_if_not_installed("shiny")

  expect_no_error(
    ddhazard_app(
      quietly = T,
      n_series = 2, sim_with = "exponential", sim_fix_options = 1,
      obs_time = 30, seed = 65848, covar_range = c(-.5, .5),
      sd_intercept = .2, sd_coef = .5, est_with_model = "exp_clip_time_w_jump",
      est_with_method = "EKF", est_fix_options = 1, LR = 1,
      order = 1, denom_term = 1, fixed_terms_method = "M_step",
      use_extra_correction = F, beta = 0, alpha = 1,
      SMA_version = "woodbury", GMA_max_rep = 25,
      GMA_NR_eps = 4, more_options = F))
})
