set.seed(12)
n_vars <- 3
sims <- test_sim_func_logit(
  n_series = 200, n_vars = n_vars, t_0 = 0, t_max = 20,
  x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
  intercept_start = -3, sds = Q_true <- c(.05, rep(.2, n_vars)))
Q_true <- diag(Q_true^2)

X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)
risk_set <-
  get_risk_obj(Y = X_Y$Y, by = 1, max_T = 20,
               id = sims$res$id, is_for_discrete_model = TRUE)

sum(sims$res$event)
xtabs(~ sims$res$tstop[sims$res$event == 1])
matplot(sims$beta, type = "l", lty = 1)

Q <- diag(1, n_vars + 1)
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
  control = list(n_max = 10))

sink("tmp.txt")
options(digits = 4)
result <- PF_EM(
  Surv(tstart, tstop, event) ~ . - id,
  data = sims$res,
  max_T = 20,
  by = 1,
  id = sims$res$id,
  Q_0 = Q_0,
  Q = Q,
  control = list(N_fw_n_bw = 1e3, N_smooth = 1e3, N_first = 2e3,
                 n_threads = 7,
                 smoother = "Brier_O_N_square"
                 ),
  trace = 2)
sink()

norm(ddfit$Q - Q_true)
norm(result$Q - Q_true)

diag(ddfit$Q - Q_true)
diag(result$Q - Q_true)

var(ddfit$state_vecs[1, ] - a_0)
var(result$a_0 - a_0)


b_idx <- 1:4 #1:5
matplot(0:20, sims$betas[, b_idx], lty = 1, type = "l", xlim = c(0, 21),
        ylim = range(ddfit$state_vecs[, b_idx], sims$betas[, b_idx]))
plot_out <- plot(result, add = TRUE, qlvls = .5, lty = 2, cov_index = b_idx)
plot(result, add = TRUE, qlvls = c(),
     type = "forward_clouds", lty = 3, cov_index = b_idx)
plot(result, add = TRUE, qlvls = c(),
     type = "backward_clouds", lty = 3, cov_index = b_idx)
matplot(0:20, ddfit$state_vecs[, b_idx], lty = 1, col = "brown", type = "l", add = TRUE)

sum((plot_out$mean - sims$betas[-1, ])^2)
sum((ddfit$state_vecs[-1, ] - sims$betas[-1, ])^2)



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









# # #TODO: clean up
# PF_effective_sample_size <- asNamespace("dynamichazard")$PF_effective_sample_size
# plot(result, type = "smoothed_clouds")
# plot(result, type = "backward_clouds", qlvls = c(), lty = 2, add = TRUE)
# plot(result, type = "forward_clouds", qlvls = c(), lty = 3, add = TRUE)
# abline(h = sims$betas[1, ])
# (tmp <- PF_effective_sample_size(result))
# (tmp2 <- PF_effective_sample_size(read_to_test("local_tests/AUX_normal_approx_w_particles")))
# for(i in seq_along(tmp)){
#   tmp[[i]] <- (tmp[[i]] - tmp2[[i]]) / tmp2[[i]]
# }
# tmp
# ddfit <- ddhazard(
#   Surv(tstart, tstop, event) ~ . - id,
#   data = sims$res,
#   max_T = 10,
#   by = 1,
#   id = sims$res$id,
#   Q_0 = diag(1, 3),
#   Q = diag(1e-1, 3),
#   a_0 = sims$betas[1, ],
#   control = list(NR_eps = 1e-5))
#
# sapply(result, function(x){
#   ws <- lapply(x, "[[", "weights")
#   sapply(ws, function(z) 1 / sum(z^2))
# })
#
# matplot(0:10, sims$betas, lty = 1, type = "l", ylim = c(-5, 5), xlim = c(0, 11))
# for(i in 1:3){
#   state_est <- t(sapply(result[[i]], function(row){
#     colSums(t(row$states) * drop(row$weights))
#   }))
#
#   idx <- switch(i, "1" = 0:10, "2" = 1:11, "3" = 1:10)
#   matplot(idx, state_est, lty = i + 1, type = "l", add = TRUE)
#   matplot(idx, state_est, lty = i + 1, type = "p", add = TRUE, pch = 16 + i)
# }
# matplot(0:10, ddfit$state_vecs, lty = 1, col = "blue", type = "l", add = TRUE)
# # #
# sapply(lapply(result$backward_clouds, "[[", "parent_idx"),
#        function(n) sort(xtabs(~ n), decreasing = TRUE)[1:10])
# sapply(lapply(result$backward_clouds, "[[", "parent_idx"),
#        function(n) sort(xtabs(~ n), decreasing = TRUE)[1:10])
#
# sapply(lapply(result$backward_clouds, "[[", "weights"),
#        function(x) sort(x, decreasing = TRUE)[1:10])
# ord <- sapply(lapply(result$backward_clouds, "[[", "weights"),
#             function(x) order(x, decreasing = TRUE)[1:11])
# #
# result$smoothed_clouds[[10]]$states[, ord[, 10]]
# result$smoothed_clouds[[8]]$states[, 5740]
#
# result <- result$smoothed_clouds
