test_that("Implement log likelihood methods and add tests",
          expect_true(FALSE))

result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1, # Use by month intervals
  n_max = 10^4, eps = 10^-4,
  a_0 = rep(0, 2), Q_0 = diag(1, 2), # Initial value
  est_Q_0 = F,
  max_T = 45,
  id = head_neck_cancer$id, order_ = 1,
  verbose = 5
)

result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  n_max = 10^4, eps = 10^-4,
  a_0 = rep(0, 4), Q_0 = diag(1, 4),
  Q = diag(c(1e-2, 1e-2, 0, 0)),
  est_Q_0 = F,
  max_T = 45,
  id = head_neck_cancer$id, order_ = 2,
  verbose = 5
)

result = ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1, # Use by month intervals
  n_max = 10^4, eps = 10^-4,
  a_0 = rep(0, 2), Q_0 = diag(1, 2), # Initial value
  est_Q_0 = F,
  max_T = 45,
  id = head_neck_cancer$id, order_ = 1,
  verbose = 5,
  model = "poisson"
)

##############
set.seed(21280)
sims <- test_sim_func_logit(n_series = 10^4, n_vars = 20, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = -.2, re_draw = T)
sims$res <- as.data.frame(sims$res)

result = ddhazard(
  survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
  sims$res,
  by = 1,
  n_max = 10^4, eps = 10^-2,
  a_0 = rep(0, 21), Q_0 = diag(1, 21),
  est_Q_0 = F,
  max_T = 10,
  id = sims$res$id, order_ = 1,
  verbose = 1
)

set.seed(31280)
sims <- test_sim_func_poisson(n_series = 10^4, n_vars = 20, t_0 = 0, t_max = 10,
                              x_range = 1, x_mean = 0, re_draw = T,
                              intercept_start = -5, sds = c(.1, rep(.4, 20)))
sum(sims$res$event)

result = ddhazard(
  survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
  sims$res,
  by = 1,
  n_max = 10^4, eps = 10^-2,
  a_0 = rep(0, 21), Q_0 = diag(1, 21),
  est_Q_0 = F,
  max_T = 10,
  id = sims$res$id, order_ = 1,
  verbose = 1,
  model = "poisson"
)

