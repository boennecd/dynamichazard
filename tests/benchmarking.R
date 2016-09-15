library(testthat); library(survival); source("R/test_utils.R")

# Simulate series to work with
set.seed(2972)
sims <- test_sim_func_logit(n_series = 10^5, n_vars = 20, t_0 = 0, t_max = 10,
                            x_range = .1, x_mean = -.4, re_draw = T)
sims$res <- as.data.frame(sims$res)

design_mat <- get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res)
rist_sets <- get_risk_obj(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

arg_list <- list(X = design_mat$X, tstart = design_mat$Y[, 1],  tstop = design_mat$Y[, 2],
                 a_0 = rep(0, ncol(design_mat$X)),
                 Q_0 = diag(10, ncol(design_mat$X)), # something large
                 Q = diag(1, ncol(design_mat$X)), # something large
                 F_ = diag(1, ncol(design_mat$X)), # first order random walk
                 risk_obj = rist_sets,
                 eps = 10^-2, n_max = 10^4,
                 order_ = 1,
                 est_Q_0 = F)

library(microbenchmark)
microbenchmark(fit2 <- do.call(ddhazard_fit_cpp_prelim, arg_list),
               fit1 <- do.call(ddhazard_fit, arg_list),
               times = 10)

test_that("Expecting similar outcome with new and old method", {
  expect_equal(c(fit1$a_t_d_s), c(fit2$a_t_d_s))
  expect_equal(c(fit1$V_t_d_s), c(fit2$V_t_d_s))
  expect_equal(c(fit1$B_s), c(fit2$B_s))
  expect_equal(c(fit1$lag_one_cor), c(fit2$lag_one_cor))
  expect_equal(c(fit1$Q), c(fit2$Q))
  expect_equal(c(fit1$Q_0), c(fit2$Q_0))
})

######
# Poisson model
set.seed(582193)
sims <- test_sim_func_poisson(1e4, n_vars = 3, x_range = 1, x_mean = .5, beta_start = 1,
                              intercept_start = -5, sds = c(.1, rep(sqrt(.5), 3)))

sum(sims$res$event)

options("scipen"=100, "digits"=4)

design_mat <- get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res)
rist_sets <- get_risk_obj(design_mat$Y, by = 1, max_T = 10, id = sims$res$id, is_for_discrete_model = F)

tmp_file <- file("benchmarking.log")
sink(tmp_file)

fit <- ddhazard_fit_cpp_prelim(
  X = design_mat$X, tstart = design_mat$Y[, 1],  tstop = design_mat$Y[, 2],
  a_0 = c(-4, 1, 1, 1),
  Q_0 = diag(10, ncol(design_mat$X)),
  Q = diag(1, ncol(design_mat$X)),
  F_ = diag(1, ncol(design_mat$X)),
  risk_obj = rist_sets,
  eps = 1e-4, n_max = 10^4,
  order_ = 1, verbose = T,
  est_Q_0 = F,
  model = "poisson",
  M_step_formulation = "Fahrmier94"
)

sink()
close(tmp_file)

matplot(sims$betas, type = "l", lty = 2)
matplot(fit$a_t_d_s, type = "l", ylim = range(fit$a_t_d_s, sims$betas), lty = 1)
matplot(sims$betas, type = "l", lty = 2, add = T)
sum(sims$res$event)
