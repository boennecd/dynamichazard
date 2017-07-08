library(dynamichazard); library(testthat); library(survival)

test_sim_func_logit <- asNamespace("dynamichazard")$test_sim_func_logit

set.seed(599479)

sims <- test_sim_func_logit(n_series = 1e5, n_vars = 10, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(10),
                            is_fixed = 2:4,
                            intercept_start = -5, sds = c(.1, rep(1, 10)))

sum(sims$res$event)

options(ddhazard_max_threads = 1)

summary(microbenchmark::microbenchmark(
  the_test = result <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~
      ddFixed(x1) + ddFixed(x2) + ddFixed(x3) + x4 + x5 + x6 + x7 + x8 + x9 + x10 -
      id - tstart - tstop - event,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(10, 8),
    Q = diag(1e-2, 8),
    control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3,
                   save_data = F, save_risk_set = F,
                   fixed_terms_method = "M_step",
                   method = "GMA"),
    max_T = 10,
    id = sims$res$id, order = 1,
    verbose = F),
  times = 3))

matplot(sims$betas, type = "l", lty = 1, col = 1:11)
matplot(result$state_vecs, type = "l", lty = 2, add = TRUE,
        col = c(1, 5:11))
for(i in 1:3)
  lines(c(0, 11), rep(result$fixed_effects[i], 2), col = i + 1, lty = 2)









set.seed(599479)

sims <- test_sim_func_logit(n_series = 500, n_vars = 10, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(10),
                            is_fixed = 2:4,
                            intercept_start = -5, sds = c(.1, rep(1, 10)))

sum(sims$res$event)

sink("tmp.txt")
result <- ddhazard(
  formula = survival::Surv(tstart, tstop, event) ~
    ddFixed(x1) + ddFixed(x2) + ddFixed(x3) + x4 + x5 + x6 + x7 + x8 + x9 + x10 -
    id - tstart - tstop - event,
  data = sims$res,
  by = (by_ <- 1),
  Q_0 = diag(10, 8),
  Q = diag(1e-2, 8),
  control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^3,
                 save_data = F, save_risk_set = F,
                 fixed_terms_method = "M_step",
                 method = "GMA", debug = TRUE),
  max_T = 10,
  id = sims$res$id, order = 1,
  verbose = F)
sink()

result$fixed_effects
