set.seed(18)
n_vars <- 6
is_fixed <- 2:7
sims <- test_sim_func_logit(
  n_series = 1e5, n_vars = n_vars, t_0 = 0, t_max = 30,
  x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(n_vars),
  is_fixed = is_fixed,
  tstart_sampl_func =
    function(t_0, t_max)
      max(t_0, runif(1, -t_max, t_max -2)),
  intercept_start = -4, sds = Q_true <- c(.05, rep(.2, n_vars)))
Q_true <- diag(Q_true[-is_fixed]^2)

X_Y = get_design_matrix(Surv(tstart, tstop, event) ~ . - id, sims$res)
risk_set <-
  get_risk_obj(Y = X_Y$Y, by = 1, max_T = 20,
               id = sims$res$id, is_for_discrete_model = TRUE)

sum(sims$res$event)
# xtabs(~ sims$res$tstop[sims$res$event == 1])
matplot(sims$beta, type = "l", lty = 1)

Q <- diag(sqrt(.1), n_vars - length(is_fixed) + 1)
Q_0 <- diag(1, n_vars - length(is_fixed) + 1)
a_0 <- sims$betas[1, ]
frm <- eval(parse(text = paste0(
  "Surv(tstart, tstop, event) ~ x1 + ", paste0(
    "ddFixed(x", is_fixed - 1, ")", collapse = " + "))))

.t <- proc.time()
ddfit <- ddhazard(
  frm,
  data = sims$res,
  max_T = 30,
  by = 1,
  id = sims$res$id,
  Q_0 = Q_0,
  model = "logit",
  Q = Q,
  control = list(n_max = 5, eps = 1e-8, n_threads = 6))
proc.time() - .t

plot(ddfit)






library(profvis)
(p <- profvis({
  result <- PF_EM(
    frm,
    data = sims$res,
    max_T = 30,
    model = 'logit',
    by = 1,
    id = sims$res$id,
    Q_0 = Q_0,
    Q = Q,
    control = PF_control(N_fw_n_bw = 1000, N_smooth = 1e4, N_first = 2e3,
                         n_threads = 7,
                         method = "AUX_normal_approx_w_cloud_mean",
                         smoother = "Fearnhead_O_N",
                         n_max = 2),
    trace = 1)
}))

sink("tmp.txt")
options(digits = 4)
.t <- proc.time()
result <- PF_EM(
  frm,
  data = sims$res,
  max_T = 30,
  model = 'logit',
  by = 1,
  id = sims$res$id,
  Q_0 = Q_0,
  Q = Q,
  control = PF_control(N_fw_n_bw = 1000, N_smooth = 1e3, N_first = 2e3,
                       n_threads = 7,
                       method = "AUX_normal_approx_w_cloud_mean",
                       smoother = "Fearnhead_O_N",
                       n_max = 1),
  trace = 2)
proc.time() - .t
sink()

norm(ddfit$Q - Q_true)
norm(result$Q - Q_true)

diag(ddfit$Q - Q_true)
diag(result$Q - Q_true)

var(ddfit$state_vecs[1, -is_fixed] - a_0[-is_fixed])
var(result$a_0 - a_0[-is_fixed])

b_idx <- setdiff(1:ncol(sims$betas), is_fixed)
ylim <- range(range(ddfit$state_vecs[, b_idx], sims$betas[, b_idx]))
ylim[1] <- min(range(ddfit$state_vecs[, b_idx], sims$betas[, b_idx]), -4)
matplot(0:30, sims$betas[, b_idx], lty = 1, type = "l", xlim = c(0, 21),
        ylim = ylim)
plot_out <- plot(result, add = TRUE, lty = 2, cov_index = b_idx)
plot(result, add = TRUE, qlvls = c(),
     type = "forward_clouds", lty = 3, cov_index = b_idx)
plot(result, add = TRUE, qlvls = c(),
     type = "backward_clouds", lty = 3, cov_index = b_idx)
matplot(0:30, ddfit$state_vecs[, b_idx], lty = 1, col = "brown", type = "l", add = TRUE)

sum((plot_out$mean - sims$betas[-1, -is_fixed])^2)
sum((ddfit$state_vecs[-1, -is_fixed] - sims$betas[-1, -is_fixed])^2)

var(result$fixed_effects - sims$betas[1, is_fixed])
var(ddfit$fixed_effects - sims$betas[1, is_fixed])








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

sink("tmp.txt")
set.seed(30302129)
result <- do.call(PF_smooth, args)
sink()
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


options(digits = 4)

beta <- c(-3, 1, .25, 1)
n <- 1e4
q <- length(beta)
tmax <- 1

X <- matrix(runif(q*n, -.5, .5), ncol = q)
y <- pmin(rexp(n, exp(drop(X %*% beta))), tmax)
sum(y == tmax)

fp <- function(y, eta, tmax)
  drop((y < tmax) - exp(eta) * y)
fpp <- function(y, eta, tmax)
  drop(- exp(eta) * y)

(cur <- beta + rnorm(length(beta), sd = sqrt(.1)))
# cur <- beta

eta <- X %*% cur
tmp <- eta  -
  fp(y = y, eta = eta, tmax = tmax) /
  fpp(y = y, eta = eta, tmax = tmax)

tmp <- t(X) %*% (tmp * (-fpp(y = y, eta = eta, tmax = tmax)))

.cov <- solve(t(X) %*% (X * (- fpp(y = y, eta = eta, tmax = tmax))))
# .cov <- t(X) %*% (X / (- fpp(y = y, eta = eta, tmax = tmax)))
.cov %*% tmp



y <- drop(1 / (1 + exp(-X %*% beta))) > runif(n)
sum(y)

fp <- function(y, eta, tmax)
  drop((exp(eta) * (y - 1) + y) / (1 + exp(eta)))
fpp <- function(y, eta, tmax)
  drop(- exp(eta) / (exp(eta) + 1)^2)

(cur <- beta + rnorm(length(beta), sd = sqrt(.1)))
# cur <- beta

eta <- X %*% cur
tmp <- eta  -
  fp(y = y, eta = eta, tmax = tmax) /
  fpp(y = y, eta = eta, tmax = tmax)

tmp <- t(X) %*% (tmp * (-fpp(y = y, eta = eta, tmax = tmax)))

.cov <- solve(t(X) %*% (X * (- fpp(y = y, eta = eta, tmax = tmax))))
# .cov <- t(X) %*% (X / (- fpp(y = y, eta = eta, tmax = tmax)))
.cov %*% tmp



.t <- 1
plot(function(x) x - exp(x) * .t, xlim = c(-200, 10), ylim = c(-170, 0))
abline(a = 0, b = 1, lty = 2)
plot(function(x) - exp(x) * .t, xlim = c(-200, 10), add = TRUE, col = "red")
abline(h = -150, tly = 3)

eta <- 10
(out <- trunc_lp_in_exponential_dist_test(eta, at_risk_length = .t, is_event = FALSE))
if(out$did_truncate)
  abline(v = out$eta_trunc, lty = 3)

out$eta_trunc - exp(out$eta_trunc) * .t
