source("C://Users//boennecd//Dropbox//skole_backup//phd//dynamichazard//R//test_utils.R")

library(dynamichazard); library(testthat); library(survival); source("C://Users//boennecd//Dropbox//skole_backup//phd//dynamichazard//R//test_utils.R")

set.seed(599479)

sims <- test_sim_func_logit(n_series = 1e6, n_vars = 10, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -5, sds = c(.1, rep(1, 10)))

getOption("ddhazard_max_threads")

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
  verbose = F))
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



set.seed(4296745)
sims <-
  test_sim_func_logit(
    n_series = 1e6, n_vars = 1, beta_start = 1,
    intercept_start = - 5, sds = c(sqrt(.1), rep(1, 1)),
    x_range = 1, x_mean = .5)$res

bench <- microbenchmark::microbenchmark(
  non_parallel =
    with(sims, get_risk_obj(
      Y = Surv(tstart, tstop, event),
      by = 1, max_T = 10, id = id, is_for_discrete_model = T)),

  parallel = with(sims, get_risk_obj(
    Y = Surv(tstart, tstop, event),
    by = 1, max_T = 10, id = id, is_for_discrete_model = T,
    n_threads = 7, min_id_size = 1e5)),

  times = 2)

summary(bench)


#####
set.seed(4296745)
sims <-
  test_sim_func_logit(
    n_series = 1e4, n_vars = 20, beta_start = rnorm(20),
    intercept_start = - 5, sds = c(sqrt(.1), rep(.3, 20)),
    x_range = 2, x_mean = .5)$res

library(profvis)
library(dynamichazard)

p <- profvis({
  tmp_mod = static_glm(Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                       data = sims, id = sims$id, by = 1,
                       control = stats::glm.control(epsilon = Inf), family = "binomial")
})

p

summary(microbenchmark::microbenchmark(
  glm = static_glm(Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                   data = sims, id = sims$id, by = 1,
                   control = stats::glm.control(epsilon = Inf), family = "binomial"),
  speedglm = static_glm(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims, id = sims$id, by = 1, family = "binomial", speedglm = T,
    maxit = 1),
  times = 5
))


summary(microbenchmark::microbenchmark(
  EKF = try(ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims, id = sims$id, by = 1,
    Q_0 = diag(1e6, 21), Q = diag(1e-2, 21),
    control = list(method = "EKF"))),

  UKF = try(ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims, id = sims$id, by = 1,
    Q_0 = diag(1, 21), Q = diag(1e-2, 21),
    control = list(method = "UKF"))),

  SMA = try(ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims, id = sims$id, by = 1,
    Q_0 = diag(1e6, 21), Q = diag(1e-2, 21),
    control = list(method = "post_approx"))),
  times = 5
))


p <- profvis({
  dd_fit <- ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims, id = sims$id, by = 1,
    Q_0 = diag(1, 21), Q = diag(1e-2, 21),
    control = list(method = "UKF", n_max = 10))
})

p

summary(microbenchmark::microbenchmark(
  UKF =  ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims, id = sims$id, by = 1,
    Q_0 = diag(1, 21), Q = diag(1e-2, 21),
    control = list(method = "UKF", n_max = 10)),
  times = 5
))
