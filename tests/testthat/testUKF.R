set.seed(2972)
sims <- test_sim_func_logit(n_series = 5e2, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = -.5, re_draw = T, beta_start = 1,
                            intercept_start = -1, sds = c(.1, rep(1, 3)))
sum(sims$res$event)

test_that("UKF does not fail",{
  res <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                  by = 1,
                  data = sims$res,
                  a_0 = rep(0, ncol(sims$res) + 1 - 4),
                  Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                  est_Q_0 = F, method = "UKF",
                  verbose = F, kappa = 0, eps = 1e-2,
                  id = sims$res$id,
                  max_T = 10)
})

EKF_res <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                    by = 1,
                    data = sims$res,
                    a_0 = rep(0, ncol(sims$res) + 1 - 4),
                    Q_0 = diag(rep(1, ncol(sims$res) + 1 - 4)),
                    est_Q_0 = F, method = "EKF",
                    verbose = F, kappa = 0,
                    id = sims$res$id,
                    max_T = 10)

matplot(res$a_t_d_s, type = "l", lty = 2, ylim = range(res$a_t_d_s, sims$betas, EKF_res$a_t_d_s))
matplot(EKF_res$a_t_d_s, type = "p", lty = 5, add = T)
matplot(sims$betas, type = "l", lty = 1, add = T)

res$Q
EKF_res$Q


#
# res_UKF <- ddhazard(formula = survival::Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
#                     by = 1,
#                     data = sims$res,
#                     a_0 = c(-1, rep(0, ncol(sims$res) + 1 - 4 - 1)),
#                     Q_0 = diag( rep(1e1, ncol(sims$res) + 1 - 4)),
#                     Q =   diag( rep(1e-1, ncol(sims$res) + 1 - 4)),
#                     est_Q_0 = F, method = "UKF", kappa = -2)
#
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
