library(dynamichazard); library(testthat); library(survival); source("C://Users//boennecd//Dropbox//skole_backup//phd//dynamichazard//R//test_utils.R")

set.seed(599479)

sims <- test_sim_func_exp(n_series = 1e6, n_vars = 10, t_0 = 0, t_max = 10,
                          x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                          intercept_start = -5, sds = c(.1, rep(1, 10)))

suppressMessages(result_exp <- ddhazard(
  formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
  data = sims$res,
  by = (by_ <- 1),
  Q_0 = diag(10, 11),
  Q = diag(1e-2, 11),
  control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3,
                 save_data = F, save_risk_set = F),
  max_T = 10,
  id = sims$res$id, order = 1,
  verbose = F,
  model = "exponential"))

matplot(result_exp$state_vecs, lty = 2, type = "l")
matplot(sims$betas, type = "l", lty = 1, add = T)
