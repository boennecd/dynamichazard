test_that("The head_neck_cancer must be defined for the gen_kalman_filter tests",
          expect_true(exists("head_neck_cancer")))

# Find the risk sets and design matrix for head and neck cancer set
design_mat <- get_design_matrix(Surv(start, stop, event) ~ group, head_neck_cancer)
rist_sets <- get_risk_sets(design_mat$Y, by = 1, max_T = 40, id = head_neck_cancer$id)

res <- ddhazard_fit(X = design_mat$X, Y = design_mat$Y,
                    a_0 = rep(0, ncol(design_mat$X)),
                    Q_0 = diag(10, ncol(design_mat$X)), # something large
                    Q = diag(1, ncol(design_mat$X)), # something large
                    F_= diag(1, ncol(design_mat$X)), # first order random walk
                    risk_sets = rist_set,
                    eps = 10^-3,
                    order_ = 1,
                    est_Q_0 = T)








  # a_0 = rep(0, ncol(design_mat$X)),
  #                            Q_0 = diag(10, ncol(design_mat$X)), # something large
  #                            Q = diag(1, ncol(design_mat$X)), # something large
  #                            F_= diag(1, ncol(design_mat$X)), # first order random walk
  #                            risk_sets = rist_set$risk_sets,
  #                            I_len = rist_set$I_len,
  #                            d = rist_set$d,
  #                            X = design_mat$X,
  #                            start = design_mat$Y[, 1],
  #                            stop = design_mat$Y[, 2],
  #                            events = design_mat$Y[, 3],
  #                            order_ = 1)
