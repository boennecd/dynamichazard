source("C://Users//boennecd//Dropbox//skole_backup//phd//dynamichazard//R//test_utils.R")
library(profvis)

l <- profvis(test_sim_func_exp(n_series = 1e4, n_vars = 10, t_0 = 0, t_max = 10,
                               x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                               intercept_start = -5, sds = c(.1, rep(1, 10))))

l

microbenchmark::microbenchmark(
  test_sim_func_exp(n_series = 2^12, n_vars = 10, t_0 = 0, t_max = 10,
                    x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                    intercept_start = -5, sds = c(.1, rep(1, 10))),
  times = 10)

# library(dynamichazard); library(testthat); library(survival); source("C://Users//boennecd//Dropbox//skole_backup//phd//dynamichazard//R//test_utils.R")
#
# set.seed(599479)
#
# sims <- test_sim_func_exp(n_series = 1e5, n_vars = 10, t_0 = 0, t_max = 10,
#                           x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
#                           intercept_start = -5, sds = c(.1, rep(1, 10)))
#
# getOption("ddhazard_max_threads")
#
# suppressMessages(result_exp <- ddhazard(
#   formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
#   data = sims$res,
#   by = (by_ <- 1),
#   Q_0 = diag(10, 11),
#   Q = diag(1e-2, 11),
#   control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3,
#                  save_data = F, save_risk_set = F),
#   max_T = 10,
#   id = sims$res$id, order = 1,
#   verbose = F,
#   model = "exponential"))
#
# matplot(result_exp$state_vecs, lty = 2, type = "l")
# matplot(sims$betas, type = "l", lty = 1, add = T)
#
#
# set.seed(45219272)
#
# sims <- test_sim_func_logit(n_series = 5e2, n_vars = 3, t_0 = 0, t_max = 10,
#                             x_range = 1, x_mean = 0.5, re_draw = T, beta_start = 2,
#                             intercept_start = -5, sds = rep(1, 2),
#                             is_fixed = 1:2)
# sum(sims$res$event)
#
#
# fit <- ddhazard(Surv(tstart, tstop, event) ~ ddFixed(1) + ddFixed(x1) + x2 + x3,
#                 data = sims$res, id = sims$res$id, by = 1, max_T = 10,
#                 model = "logit", Q_0 = diag(10, 2), Q = diag(.01, 2))
#
# matplot(sims$betas, lty = 1, type = "l", ylim = range(sims$betas, fit$state_vecs, fit$fixed_effects))
# matplot(fit$state_vecs, lty = 2, type = "l", add = T, col = 3:4)
# abline(h = fit$fixed_effects, lty = 2, col = 1:2)
#
# fit <- ddhazard(Surv(tstart, tstop, event) ~ ddFixed(1) + ddFixed(x1) + x2 + x3,
#                 data = sims$res, id = sims$res$id, by = 1, max_T = 10,
#                 model = "logit", Q_0 = diag(10, 2), Q = diag(.01, 2),
#                 control = list(fixed_terms_method = "E_step"))
#
# matplot(fit$state_vecs, lty = 3, type = "l", add = T, col = 3:4)
# abline(h = fit$fixed_effects, lty = 3, col = 1:2)
