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



set.seed(65848L)
dat <- test_sim_func_logit(n_series = 65536, n_vars = 5,
                           t_max = 30, re_draw = T, beta_start = runif(5,
                                                                       min = -1.5, max = 1.5), intercept_start = -3.5,
                           sds = c(0.4, rep(1L, 5)), x_range = 3L, x_mean = -0.5,
                           lambda = 0.166666666666667, is_fixed = NULL,
                           tstart_sampl_func = function(t_0 = t_0, t_max = t_max) max(0,
                                                                                      runif(1, t_0 - t_max, t_max - 1 - 1e-08)))
dat$res[dat$res$event == 0 & dat$res$tstop > 30,
        "tstop"] <- 30

tmp <- file("tmp.txt")
sink(tmp)
fit <- ddhazard(formula = survival::Surv(tstart,
                                  tstop, event) ~ x1 + x2 + x3 +
           x4 + x5, data = dat$res, by = 1,
         Q_0 = diag(1, 6), Q = diag(0.1,
                                    6), max_T = 30L, id = dat$res$id,
         order = 1L, model = "logit",
         control = list(eps = 0.01,
                        method = "GMA", GMA_max_rep = 10L,
                        GMA_NR_eps = 0.1, debug = T,
                        LR = .5))
sink()
close(tmp)


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
    n_series = 2e5, n_vars = 4, beta_start = rnorm(4),
    intercept_start = - 5, sds = c(sqrt(.1), rep(.3, 4)),
    x_range = 2, x_mean = .5)$res

library(profvis)
library(dynamichazard)

(p <- profvis({
  tmp_mod = static_glm(Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                       data = sims, id = sims$id, by = 1,
                       control = stats::glm.control(epsilon = Inf), family = "binomial")
}))


frm <- Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id

tmp <- get_design_matrix(frm, sims)

(p <- profvis({
  tmp_mod = static_glm(Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
                       data = sims, id = sims$id, by = 1,
                       epsilon = Inf, family = "binomial",
                       only_coef = TRUE, mf = cbind(tmp$X, tmp$fixed_terms))
}))

(p <- profvis({
  static_glm(
    frm,
    data = sims, id = sims$id, by = 1, family = "binomial", speedglm = T,
    maxit = 1, only_coef = TRUE, mf = cbind(tmp$X, tmp$fixed_terms))
}))


summary(microbenchmark::microbenchmark(
   glm = (f1 <- static_glm(frm,
                   data = sims, id = sims$id, by = 1,
                   control = stats::glm.control(epsilon = Inf), family = "binomial")),

  speedglm = (f2 <- static_glm(
    frm,
    data = sims, id = sims$id, by = 1, family = "binomial", speedglm = T,
    maxit = 1)),

  glm_coef = (f3 <-
                static_glm(frm,
                           data = sims, id = sims$id, by = 1,
                           epsilon = Inf, family = "binomial",
                           only_coef = TRUE, mf = cbind(tmp$X, tmp$fixed_terms))),

  speedglm_coef = (f4 <- static_glm(
    frm,
    data = sims, id = sims$id, by = 1, family = "binomial", speedglm = T,
    maxit = 1, only_coef = TRUE, mf = cbind(tmp$X, tmp$fixed_terms))),


  times = 10
))

rbind(f1$coefficients, f2$coefficients, f3, f4)


summary(microbenchmark::microbenchmark(
  EKF = try(ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims, id = sims$id, by = 1,
    Q_0 = diag(1, 5), Q = diag(1e-2, 5),
    control = list(method = "GMA"))),
  times = 25
))


set.seed(4296745)
sims <-
  test_sim_func_logit(
    n_series = 1e5, n_vars = 4, beta_start = rnorm(4),
    intercept_start = - 3, sds = c(sqrt(.1), rep(.5, 4)),
    x_range = 2, x_mean = 0)
sum(sims$res$event)

options(ddhazard_use_speedglm = T)

p <- profvis({
  dd_fit <- ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims$res, id = sims$res$id, by = 1,
    Q_0 = diag(1, 5), Q = diag(1e-2, 5),
    control = list(method = "EKF", n_max = 10,
                   NR_eps = .1))
})

p

summary(microbenchmark::microbenchmark(
  EKF = ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims$res, id = sims$res$id, by = 1,
    Q_0 = diag(1, 5), Q = diag(1e-2, 5),
    control = list(method = "EKF", n_max = 10,
                   NR_eps = .1, EKF_batch_size = 2500L)),

  EKF_smaller = ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims$res, id = sims$res$id, by = 1,
    Q_0 = diag(1, 5), Q = diag(1e-2, 5),
    control = list(method = "EKF", n_max = 10,
                   NR_eps = .1, EKF_batch_size = 1000L)),

  EKF_smallest = ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims$res, id = sims$res$id, by = 1,
    Q_0 = diag(1, 5), Q = diag(1e-2, 5),
    control = list(method = "EKF", n_max = 10,
                   NR_eps = .1, EKF_batch_size = 250L)),

  GMA = ddhazard(
    Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
    data = sims$res, id = sims$res$id, by = 1,
    Q_0 = diag(1, 5), Q = diag(1e-2, 5),
    control = list(method = "GMA", n_max = 10,
                   GMA_NR_eps = .1)),

  times = 20
))
