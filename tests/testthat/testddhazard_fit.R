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

test_that("Implement tests",
          expect_true(FALSE))
