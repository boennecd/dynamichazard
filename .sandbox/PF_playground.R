PF_smooth <- asNamespace("dynamichazard")$PF_smooth

n_vars <- 3
set.seed(78095325)
sims <- test_sim_func_logit(
  n_series = 1e4, n_vars = n_vars, t_0 = 0, t_max = 10,
  x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
  intercept_start = -3, sds = c(.1, rep(.5, n_vars)))

X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)
risk_set <-
  get_risk_obj(Y = X_Y$Y, by = 1, max_T = 10,
               id = sims$res$id, is_for_discrete_model = TRUE)

sum(sims$res$event)
matplot(sims$beta, type = "l", lty = 1)

Q <- diag(c(.1, rep(.5, n_vars)))
Q_0 <- diag(1, n_vars + 1)
a_0 <- sims$betas[1, ]

args <- list(
  n_fixed_terms_in_state_vec = 0,
  X = t(X_Y$X),
  fixed_terms = t(X_Y$fixed_terms),
  tstart = X_Y$Y[1, ],
  tstop = X_Y$Y[2, ],
  Q_0 = Q_0,

  Q = Q,
  a_0 = a_0,


  Q_tilde = diag(0, n_vars + 1),
  risk_obj = risk_set,
  F = diag(1, n_vars + 1),
  n_max = 10,
  order = 1,
  n_threads = 7,
  N_fw_n_bw = 1e3,
  N_smooth = 1e3,
  N_first = 1e4,
  forward_backward_ESS_threshold = NULL,
  debug = 2, # TODO: remove
  method = "PF")

sink("tmp.txt") # TODO: remove
set.seed(30302129)
old_seed <- .Random.seed
args$method <-  "AUX_normal_approx_w_particles"  # "AUX_temp"
result <- do.call(PF_smooth, args)
sink() # TODO: remove

sapply(result, function(x){
  ws <- lapply(x, "[[", "weights")
  sapply(ws, function(z) 1 / sum(z^2))
})




ddfit <- ddhazard(
  Surv(tstart, tstop, event) ~ . - id,
  data = sims$res,
  max_T = 10,
  by = 1,
  id = sims$res$id,
  Q_0 = Q_0,
  Q = Q,
  a_0 = sims$betas[1, ])


b_idx <- 1:4 #1:5
matplot(0:10, sims$betas[, b_idx], lty = 1, type = "l", ylim = c(-5, 5), xlim = c(0, 11))
for(i in 2:3){
  state_est <- t(sapply(result[[i]], function(row){
    colSums(t(row$states) * drop(row$weights))
  }))

  idx <- switch(i, "1" = 0:10, "2" = 1:11, "3" = 1:10)
  matplot(idx, state_est[, b_idx], lty = i + 1, type = "l", add = TRUE)
  matplot(idx, state_est[, b_idx], lty = i + 1, type = "p", add = TRUE, pch = 16 + i)
}
matplot(0:10, ddfit$state_vecs[, b_idx], lty = 1, col = "brown", type = "l", add = TRUE)

sum((state_est - sims$betas[-1, ])^2)
sum((ddfit$state_vecs[-1, ] - sims$betas[-1, ])^2)

options(digits = 3)
sapply(lapply(result$smoothed_clouds, "[[", "weights"),
       function(x) sort(x, decreasing = TRUE)[1:10])
ord <- sapply(lapply(result$smoothed_clouds, "[[", "weights"),
            function(x) order(x, decreasing = TRUE)[1:100])

t <- 5
result$smoothed_clouds[[t]]$states[, ord[1:20, t]]
sims$betas[t + 1, ]

colSums((
  result$smoothed_clouds[[t]]$states[, ord[1:20, t]] -
    sims$betas[t + 1, ])^2)
