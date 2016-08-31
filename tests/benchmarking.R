library(testthat); library(survival); library(benssurvutils); source("R/test_utils.R")

# Simulate series to work with
set.seed(2972)
sims <- test_sim_func_logit(n_series = 10^5, n_vars = 20, t_0 = 0, t_max = 10,
                            x_range = .1, x_mean = -.4, re_draw = T)
sims$res <- as.data.frame(sims$res)

design_mat <- benssurvutils::get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res)
rist_sets <- benssurvutils::get_risk_sets(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

design_mat$Y[, 2] <- rist_sets$stop_new
design_mat$Y[, 3] <- rist_sets$new_events_flags

arg_list <- list(X = design_mat$X, tstart = design_mat$Y[, 1],  tstop = design_mat$Y[, 2], events = design_mat$Y[, 3],
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
               times = 1)

test_that("Expecting similar outcome with new and old method", {
  expect_equal(c(fit1$a_t_d_s), c(fit2$a_t_d_s))
  expect_equal(c(fit1$V_t_d_s), c(fit2$V_t_d_s))
  expect_equal(c(fit1$B_s), c(fit2$B_s))
  expect_equal(c(fit1$lag_one_cor), c(fit2$lag_one_cor))
  expect_equal(c(fit1$Q), c(fit2$Q))
  expect_equal(c(fit1$Q_0), c(fit2$Q_0))
})
