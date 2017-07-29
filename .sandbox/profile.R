library(dynamichazard);
library(survival);
test_sim_func_logit <- asNamespace("dynamichazard")$test_sim_func_logit

set.seed(9997)
sims <- test_sim_func_logit(n_series = 1e5, n_vars = 10, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = (1:10 - 5) / 2.5,
                            intercept_start = -5, sds = c(.1, rep(1, 10)))

microbenchmark::microbenchmark(
  exp = fit <- ddhazard(
    formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    data = sims$res,
    by = (by_ <- 1),
    Q_0 = diag(1, 11),
    Q = diag(1e-1, 11),
    control = list(est_Q_0 = F, eps = 10^-3,
                   method = "EKF", NR_eps = 1e-5),
    max_T = 10,
    id = sims$res$id, order = 1),
  times = 5)

# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# exp 4.546921 4.619881 4.713849 4.720406 4.784067 4.897973     5
