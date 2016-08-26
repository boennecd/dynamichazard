# library(survival); library(benssurvutils); library(dynamichazard); source("R/test_utils.R")

# Simulate series to work with
set.seed(2972)
sims <- test_sim_func_logit(n_series = 10^4, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = .1, x_mean = -.4, re_draw = T)
sims$res <- as.data.frame(sims$res)

design_mat <- benssurvutils::get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res)
rist_sets <- benssurvutils::get_risk_sets(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

design_mat$Y[, 2] <- rist_sets$stop_new
design_mat$Y[, 3] <- rist_sets$new_events_flags

# tmp_file <- file("testddhazard_fit_R_vs_cpp_1.log")
# sink(tmp_file)
tmp_time <- Sys.time()
res <- ddhazard_fit(X = design_mat$X, Y = design_mat$Y,
                    a_0 = rep(0, ncol(design_mat$X)),
                    Q_0 = diag(10, ncol(design_mat$X)), # something large
                    Q = diag(1, ncol(design_mat$X)), # something large
                    F_= diag(1, ncol(design_mat$X)), # first order random walk
                    risk_sets = rist_sets,
                    eps = 10^-4, n_max = 10^4,
                    order_ = 1,
                    est_Q_0 = F)
Sys.time() - tmp_time
# sink()
# close(tmp_file)

# tmp_file <- file("testddhazard_fit_R_vs_cpp_2.log")
# sink(tmp_file)
tmp_time <- Sys.time()
res_new <- ddhazard_fit_cpp_prelim(
  X = design_mat$X,
  tstart = design_mat$Y[, 1],  tstop = design_mat$Y[, 2], events = design_mat$Y[, 3],
  a_0 = rep(0, ncol(design_mat$X)),
  Q_0 = diag(10, ncol(design_mat$X)), # something large
  Q = diag(1, ncol(design_mat$X)), # something large
  F = diag(1, ncol(design_mat$X)), # first order random walk
  risk_obj = rist_sets,
  eps = 10^-4, n_max = 10^4,
  order = 1,
  est_Q_0 = F)
Sys.time() - tmp_time
# sink()
# close(tmp_file)

test_that("Expecting similar outcome with new and old method", {
  expect_equal(c(res_new$a_t_d_s), c(res$a_t_d_s))
})

test_that("Check other indicies",
          expect_true(F))

test_that("Check other order, with and without estimating Q_0",
          expect_true(F))

# par(mfcol = c(2, 2))
# plot(res$a_t_d_s[, 1], type = "l", ylim = range(res$a_t_d_s[, 1], res_new$a_t_d_s[, 1]))
# lines(seq_len(nrow(sims$betas) + 1), res_new$a_t_d_s[, 1], col = "blue")
# for(i in 1:3){
#   plot(res$a_t_d_s[, i + 1], type = "l", ylim = range(sims$betas[, i ], res$a_t_d_s[, i + 1], res_new$a_t_d_s[, i + 1]))
#   lines(x = seq_len(nrow(sims$betas)), sims$betas[, i], col = "red")
#   lines(x = seq_len(nrow(sims$betas) + 1), res_new$a_t_d_s[, i + 1], col = "blue")
# }
