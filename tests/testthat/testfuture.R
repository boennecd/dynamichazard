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
  set.seed(10)
  ws <- sample.int(nrow(head_neck_cance), 60)
  ws <- sapply(1:nrow(head_neck_cancer), function(x) sum(ws == x))

  expect_silent(
    ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 3,
      a_0 = rep(0, 2), Q_0 = diag(1, 2),
      control = list(method = meth),
      max_T = 27, order = 1,
      weights = ws))
})

test_that("Making large design mat and using weights yield the same",{
  set.seed(9191)
  tmp <- sample.int(nrow(head_neck_cancer), 100, replace = T)
  dum_design <- head_neck_cancer[tmp, ]

  ws <- sapply(1:nrow(head_neck_cancer), function(x) sum(tmp == x))

  meth <- "UKF"
  for(m in c("logit")){
    f1 <- ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = dum_design,
      by = 5, model = m,
      a_0 = c(-2,0), Q_0 = diag(1, 2), Q = diag(1e-2, 2),
      control = list(method = meth, ridge_eps = 1e-3),
      max_T = 25, order = 1)

    suppressMessages(f2 <- ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 5, model = m,
      a_0 = c(-2, 0), Q_0 = diag(1, 2), Q = diag(1e-2, 2),
      control = list(method = meth, ridge_eps = 1e-3),
      max_T = 25, order = 1,
      weights = ws))

    info <- paste("m =", m)
    expect_equal(f1$state_vecs, f2$state_vecs, info = info, tolerance = 1e-5)
    # expect_equal(f1$state_vars, f2$state_vars, info = info, tolerance = 1e-5)
  }
})


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
