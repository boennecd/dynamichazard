set.seed(2972)
sims <- test_sim_func_logit(n_series = 5e2, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = -.5, re_draw = T, beta_start = 1,
                            intercept_start = -1, sds = c(.1, rep(1, 3)))
sum(sims$res$event)

test_that("UKF does not fail and both methods give the same",{
  res_new <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                  by = 1,
                  data = sims$res,
                  a_0 = rep(0, ncol(sims$res) + 1 - 4),
                  Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                  est_Q_0 = F, method = "UKF",
                  verbose = F, kappa = 0, eps = 1e-2,
                  id = sims$res$id,
                  max_T = 10)


  res_old <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                      by = 1,
                      data = sims$res,
                      a_0 = rep(0, ncol(sims$res) + 1 - 4),
                      Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                      est_Q_0 = F, method = "UKF_org",
                      verbose = F, kappa = 0, eps = 1e-2,
                      id = sims$res$id,
                      max_T = 10)

  expect_equal(res_new$a_t_d_s, res_old$a_t_d_s)
  expect_equal(res_new$V_t_d_s, res_old$V_t_d_s)
})

# res_test <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
#                      by = 1,
#                      data = sims$res,
#                      a_0 = res_new$a_t_d_s[1, ],
#                      Q_0 = res_new$Q,
#                      Q = res_new$Q,
#                      est_Q_0 = F, method = "UKF", kappa = NULL)

# res_EKF <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
#                     by = 1,
#                     data = sims$res,
#                     a_0 = rep(0, ncol(sims$res) + 1 -4),
#                     Q_0 = diag(rep(1, ncol(sims$res) + 1 -4)),
#                     est_Q_0 = F)
#
# matplot(res_UKF$a_t_d_s, ylim = range(res_UKF$a_t_d_s, res_EKF$a_t_d_s),
#         lty = 1, type = "l")
# matplot(res_EKF$a_t_d_s, lty = 2, type = "l", add = T)

test_that("Implement tests for UKF", {
  expect_true(F)
})
