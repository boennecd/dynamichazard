PF_smooth <- asNamespace("dynamichazard")$PF_smooth

set.seed(5)
n_vars <- 3
sims <- test_sim_func_logit(
  n_series = 200, n_vars = n_vars, t_0 = 0, t_max = 20,
  x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
  intercept_start = -3, sds = c(.05, rep(.2, n_vars)))

X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)
risk_set <-
  get_risk_obj(Y = X_Y$Y, by = 1, max_T = 20,
               id = sims$res$id, is_for_discrete_model = TRUE)

sum(sims$res$event)
xtabs(~ sims$res$tstop[sims$res$event == 1])
matplot(sims$beta, type = "l", lty = 1)

Q <- diag(.5, n_vars + 1)
Q_0 <- diag(.1, n_vars + 1)
a_0 <- sims$betas[1, ]

ddfit <- ddhazard(
  Surv(tstart, tstop, event) ~ . - id,
  data = sims$res,
  max_T = 20,
  by = 1,
  id = sims$res$id,
  Q_0 = Q_0,
  Q = Q,
  a_0 = sims$betas[1, ])

sink("tmp.txt")
tmp <- PF_EM(
  n_fixed_terms_in_state_vec = 0,
  X = t(X_Y$X),
  fixed_terms = t(X_Y$fixed_terms),
  tstart = X_Y$Y[1, ],
  tstop = X_Y$Y[2, ],
  Q_0 = Q_0,

  Q = Q,
  a_0 = a_0,
  risk_obj = risk_set,
  n_max = 200,
  order = 1,
  n_threads = 7,
  N_fw_n_bw = 1e3,
  N_smooth = 1e3,
  N_first = 1e4,
  forward_backward_ESS_threshold = NULL,
  trace = 3, # TODO: remove
  method = "AUX_normal_approx_w_particles",
  eps = 1e-2)
sink()



print(loglikes <- c(loglikes, logLik(result)), digits = 6)
PF_effective_sample_size(result)

b_idx <- 1:4 #1:5
matplot(0:20, sims$betas[, b_idx], lty = 1, type = "l", ylim = c(-5, 5), xlim = c(0, 21),
        main = .i)
for(i in 3){
  state_est <- t(sapply(result[[i]], function(row){
    colSums(t(row$states) * drop(row$weights))
  }))

  idx <- switch(i, "1" = 0:20, "2" = 1:21, "3" = 1:20)
  matplot(idx, state_est[, b_idx], lty = i + 1, type = "l", add = TRUE, lwd = if(3 == i) 1.5 else .5)
  matplot(idx, state_est[, b_idx], lty = i + 1, type = "p", add = TRUE, pch = 16 + i,
          cex = if(3 == i) 1 else .75)
}
matplot(0:20, ddfit$state_vecs[, b_idx], lty = 1, col = "brown", type = "l", add = TRUE)

sts <- compute_summary_stats(result)

a_0 <- drop(sts[[1]]$E_xs)
Q <- matrix(0., length(a_0), length(a_0))
for(i in 1:length(sts))
  Q <- Q + sts[[i]]$E_x_less_x_less_one_outers
Q <- Q / length(sts)

args$a_0 <- a_0
args$Q <- Q

a_0  - ddfit$state_vecs[1, ]
diag(Q) - diag(ddfit$Q)

sink("tmp.txt") # TODO: remove
set.seed(30302129)
result <- do.call(PF_smooth, args)
sink() # TODO: remove

.i <- .i + 1







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
