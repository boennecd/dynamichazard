test_that("The head_neck_cancer must be defined for the ddhazard_fit tests",
          expect_true(exists("head_neck_cancer")))

# Find the risk sets and design matrix for head and neck cancer set
design_mat <- benssurvutils::get_design_matrix(Surv(start, stop, event) ~ group, head_neck_cancer)
rist_sets <- benssurvutils::get_risk_sets(design_mat$Y, by = 1, max_T = 45, id = head_neck_cancer$id)

res <- ddhazard_fit(X = design_mat$X, Y = design_mat$Y,
                    a_0 = rep(0, ncol(design_mat$X)),
                    Q_0 = diag(1, ncol(design_mat$X)), # something large
                    Q = diag(1, ncol(design_mat$X)), # something large
                    F_= diag(1, ncol(design_mat$X)), # first order random walk
                    risk_sets = rist_sets,
                    eps = 10^-4,
                    order_ = 1,
                    est_Q_0 = T)

# Simulate series to work with
set.seed(2972)
sims <- test_sim_func_logit(n_series = 10^3, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = -.1, re_draw = T)
sims$res <- as.data.frame(sims$res)

design_mat <- benssurvutils::get_design_matrix(Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res)
design_mat$X <- design_mat$X[, -1] # remove intercept
rist_sets <- benssurvutils::get_risk_sets(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

res <- ddhazard_fit(X = design_mat$X, Y = design_mat$Y,
                    a_0 = rep(0, ncol(design_mat$X)),
                    Q_0 = diag(1, ncol(design_mat$X)), # something large
                    Q = diag(1, ncol(design_mat$X)), # something large
                    F_= diag(1, ncol(design_mat$X)), # first order random walk
                    risk_sets = rist_sets,
                    eps = 10^-4, n_max = 10^4,
                    order_ = 1,
                    est_Q_0 = F)

res

plot(res$a_t_d_s[, 1], type = "l")
plot(sims$betas[, 1], type = "l")

plot(res$a_t_d_s[, 2], type = "l")
plot(sims$betas[, 2], type = "l")

plot(res$a_t_d_s[, 3], type = "l")
plot(sims$betas[, 3], type = "l")
