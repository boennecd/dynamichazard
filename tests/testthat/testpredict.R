# library(testthat); library(survival); library(parallel); source("R/test_utils.R")

result = ddhazard(
  formula = Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2), Q_0 = diag(.1, 2),
  save_risk_set = T
)

for(use_parallel in c(T, F)){
  # Test that we get same estimate in one period estimates
  test_that(paste0("Testing one period in sample with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(1, 60)),
                       use_parallel = use_parallel)
    tmp_a_t_d_s <- rbind(result$a_t_d_s[-1, ], result$a_t_d_s[60, ])
    expect_equal(predict_$fits, c(result$hazard_func(tmp_a_t_d_s %*% c(1, 1))),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(0, 60)),
                       use_parallel = use_parallel)
    expect_equal(predict_$fits, c(result$hazard_func(tmp_a_t_d_s %*% c(1, 0))),
                 use.names = F, check.attributes = F)
  })

  # Check with two period estimates
  test_that(paste0("Testing two period in sample with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 2*(0:29), stop = 2*(1:30), group = rep(1, 30)),
                       use_parallel = use_parallel)
    tmp_a_t_d_s <- rbind(result$a_t_d_s[-1, ], result$a_t_d_s[60, ])
    fac1 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30) - 1, ] %*% c(1, 1)))
    fac2 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30), ] %*% c(1, 1)))
    expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 2*(0:29), stop = 2*(1:30), group = rep(0, 30)),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30) - 1, ] %*% c(1, 0)))
    fac2 = c(result$hazard_func(tmp_a_t_d_s[2*(1:30), ] %*% c(1, 0)))
    expect_equal(predict_$fits, 1 - (1 - fac1) * (1 - fac2),
                 use.names = F, check.attributes = F)
  })

  # Check forcasting
  test_that(paste0("Testing forcasting with use_parallel = ", use_parallel),{
    predict_ = predict(result, new_data = data.frame(start = 59, stop = 69, group = 1),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(result$a_t_d_s[60, ] %*% c(1, 1)))
    expect_equal(predict_$fits, (1 - (1 - fac1)^10),
                 use.names = F, check.attributes = F)

    predict_ = predict(result, new_data = data.frame(start = 59, stop = 69, group = 0),
                       use_parallel = use_parallel)
    fac1 = c(result$hazard_func(result$a_t_d_s[60, ] %*% c(1, 0)))
    expect_equal(predict_$fits, (1 - (1 - fac1)^10),
                 use.names = F, check.attributes = F)
  })
}

# Check the terms
test_that("Testing term prediction",{
  for(g in 0:1){
    predict_ = predict(result, new_data = data.frame(start = 0:58, stop = 1:59, group = g))
    predict_terms = predict(result, new_data = data.frame(group = g), type = "term")
    expect_equal(dim(predict_terms$terms), c(60, 1, 2))
    predict_terms$terms <- predict_terms$terms[-1, , , drop = F]

    predict_terms_fit = exp(c(apply(predict_terms$terms, 1:2, sum))) /
      (exp(c(apply(predict_terms$terms, 1:2, sum))) + 1)
    expect_equal(unname(predict_$fits), predict_terms_fit, check.attributes = F)
  }
})

test_that("Terms prediction when parsing two observations",{
  predict_terms = predict(result, new_data = data.frame(group = 0:1), type = "term")
  expect_equal(dim(predict_terms$terms), c(60, 2, 2))
})

# Check for second order
result = ddhazard(
  formula = Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  a_0 = rep(0, 2 * 2), Q_0 = diag(10, 2 * 2),
  Q = diag(c(1.0e-4, 1.0e-4, 0, 0)), est_Q_0 = F,
  n_max = 1e3,
  save_risk_set = T, order_ = 2
)

test_that("Calls with second order models do not throw errors", {
  for(g in c(0, 1))
    for(use_parallel in c(T, F))
      expect_no_error(bquote(
        predict(result, new_data = data.frame(start = 0:59, stop = 1:60, group = rep(.(g), 60)),
                use_parallel = .(use_parallel))))
})

dum <- structure(list(model = "poisson"), "class" = class(result))

test_that("predict functions throws error when model is poisson",{
  expect_error(predict(dum))
})

warning("Implement test for second order random walk and predict")

######
# Poisson model
set.seed(39795)
sims <- test_sim_func_poisson(1e4, n_vars = 3, x_range = 1, x_mean = .5, beta_start = 1,
                              intercept_start = -5, sds = c(.2, rep(.5, 3)))

design_mat <- get_design_matrix(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res)
risk_obj <- get_risk_obj(design_mat$Y, by = 1, max_T = 10, id = sims$res$id)

fit <- ddhazard(
  survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3, sims$res,
  a_0 = rep(0, ncol(design_mat$X)),
  Q_0 = diag(10, ncol(design_mat$X)),
  Q = diag(1, ncol(design_mat$X)),
  eps = 1e-4, n_max = 10^4,
  order_ = 1, by = 1,
  est_Q_0 = F,
  model = "poisson"
)

matplot(fit$a_t_d_s, type = "l", ylim = range(fit$a_t_d_s, sims$betas), lty = 1)
matplot(sims$betas, type = "l", lty = 2, add = T)
sum(sims$res$event)

preds <- predict(fit, new_data = sims$res, tstart = "tstart", tstop = "tstop")

test_that("poisson predicting have equal number of rows", {
  expect_equal(length(preds$fits), nrow(sims$res))
})

test_that("random sample rows have correct predicted chance of failure given estimates in Poisson model",{
  r_rows <- sample.int(nrow(sims$res), min(nrow(sims$res), 1e3), replace = F)
  l_parem <- rbind(fit$a_t_d_s[-1, ], fit$a_t_d_s[nrow(fit$a_t_d_s), ])
  l_time <- c(fit$times, max(sims$res$tstop))

  for(r in r_rows){
    x <- c(1, unlist(sims$res[r, -(1:4)]))
    tstart <- sims$res$tstart[r]
    tstop <- sims$res$tstop[r]

    istart <- max(which(l_time <= tstart))
    istop <- min(which(l_time[-1] >= tstop))

    p_surv <- 1
    for(i in istart:istop){
      p_surv = p_surv * (
        1 - fit$hazard_func(eta = (l_parem[i, ] %*% x)[1,1],
                                    tstart = max(l_time[i], tstart),
                                    tstop = min(l_time[i + 1], tstop)))
    }

    if(is.character(all.equal(1 - p_surv, unname(preds$fits[r]))))
      expect_true(FALSE, info = paste0(c("Predict failed for row with start and stop =", tstart, tstop), collapse = " "))
  }
})
