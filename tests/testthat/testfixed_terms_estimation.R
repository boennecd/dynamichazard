library(biglm)

if(interactive()){
  library(survival); library(dynamichazard); library(testthat)
  source("C:/Users/boennecd/Dropbox/skole_backup/phd/dynamichazard/R/test_utils.R")
}

set.seed(548237)
sims <- test_sim_func_logit(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -4, sds = c(.1, rep(1, 3)))
# sum(sims$res$event)

test_that("Only fixed effects yields same results as bigglm with logit model", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  suppressWarnings(
    res1 <- ddhazard(form, data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10,
                     control = list(eps_fixed_parems = 1e-3, fixed_effect_chunk_size = 1e3, max_it_fixed_parems = 10)))

  tmp_design <- get_survival_case_weigths_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10)

  suppressWarnings(res2 <- bigglm(update(form, Y ~ .), data = tmp_design, family = binomial(), chunksize = 1e3))

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects)))
})


# form <- formula(survival::Surv(tstart, tstop, event) ~
#                   -1 + ddFixed(rep(1, length(x1))) + x1 + x2 + ddFixed(x3))
#
# res1 <- ddhazard(form, data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10)
#
# matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
# matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
# abline(h = res1$fixed_effects, col = c(1,4))
#
#
#
#
#
#
#
# set.seed(548237)
# sims <- test_sim_func_logit(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
#                             x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
#                             intercept_start = -4, sds = c(.1, rep(1, 3)),
#                             is_fixed = c(1, 4))
# sum(sims$res$event)
#
# res1 <- ddhazard(form, data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10)
#
# matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
# matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
# abline(h = res1$fixed_effects, col = c(1,4))

if(!interactive()) test_that("Make the above into a test", expect_true(F))

set.seed(312237)
sims <- test_sim_func_exp(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
                          x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                          intercept_start = -4, sds = c(.1, rep(1, 3)))
# sum(sims$res$event)

test_that("Only fixed effects yields same results as bigglm with exponential model", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  suppressWarnings(res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
                                    control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3)))

  tmp_design <- get_survival_case_weigths_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10, is_for_discrete_model = F)

  suppressWarnings(res2 <- bigglm(update(form, Y ~ . + offset(log(pmin(tstop, t) - pmax(tstart, t - 1)))), data = tmp_design, family = poisson(), chunksize = 1e3))

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects)))
})


test_that("Changing fixed effect control parems changes the result", {
  arg_list <- list(
    formula(survival::Surv(tstart, tstop, event) ~
              -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3)),
    data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
    control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3))

  suppressWarnings(res1 <- do.call(ddhazard, arg_list))

  # Should not make a difference
  arg_list_tmp <- arg_list
  arg_list_tmp$control$fixed_effect_chunk_size <- 1e4
  suppressWarnings(res2 <- do.call(ddhazard, arg_list_tmp))
  expect_equal(res2$fixed_effects, res1$fixed_effects)

  # Should make a difference
  arg_list_tmp <- arg_list
  arg_list_tmp$control$eps_fixed_parems <- 1
  suppressWarnings(res2 <- do.call(ddhazard, arg_list_tmp))
  expect_true(!all(res2$fixed_effects == res1$fixed_effects))

  # Should make a difference
  arg_list_tmp <- arg_list
  arg_list_tmp$control$max_it_fixed_parems <- 5
  suppressWarnings(res2 <- do.call(ddhazard, arg_list_tmp))
  expect_true(!all(res2$fixed_effects == res1$fixed_effects))
})

#
#
#
# form <- formula(survival::Surv(tstart, tstop, event) ~
#                   -1 + ddFixed(rep(1, length(x1))) + x1 + x2 + x3)
#
# res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
#                  control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3, LR = .5), Q_0 = diag(rep(10, 3)))
#
# matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
# matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
# abline(h = res1$fixed_effects, col = c(1,4))
#
#
#
#
#
# form <- formula(survival::Surv(tstart, tstop, event) ~
#                    + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))
#
# res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10)
#
# matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
# matplot(res1$state_vecs, add = T, col = 1, type = "l", lty = 1)
# abline(h = res1$fixed_effects, col = 2:4)
#
#
#
#
#
#
# set.seed(312237)
# sims <- test_sim_func_exp(n_series = 1e5, n_vars = 3, t_0 = 0, t_max = 10,
#                           x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
#                           intercept_start = -4, sds = c(.1, rep(1, 3)), is_fixed = 2:3)
#
# form <- formula(survival::Surv(tstart, tstop, event) ~ 1 + ddFixed(x1) + ddFixed(x2) + x3)
#
# res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10)
#
# matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
# matplot(res1$state_vecs, add = T, col = c(1,4), type = "l", lty = 1)
# abline(h = res1$fixed_effects, col = 2:3)



if(!interactive()) test_that("Make the above into a test", expect_true(F))

test_that("UKF with fixed effects works", {
  set.seed(2231412)
  sims <- test_sim_func_exp(n_series = 1e4, n_vars = 3, t_0 = 0, t_max = 10,
                            x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                            intercept_start = -4, sds = c(.1, rep(1, 3)))

  fit <- ddhazard(formula(survival::Surv(tstart, tstop, event) ~
                            -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + x2 + x3),
                  data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10,
                  control = list(method = "UKF", fixed_parems_start = rep(0, 2)))


  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(fit$state_vecs, type = "l", lty = 2, col = 3:4, add = T)
  # abline(h = fit$fixed_effects, col = 1:2)
  # get_expect_equal(fit, file = "tmp.txt")

  expect_equal(names(fit$fixed_effects), c("rep(1, length(x1))", "x1"))

  expect_equal(c(fit$state_vecs),
               c( 0.007388578654419261, 0.007389254813363675, -0.008510186512430018, 0.275876744652119421, 0.650608217501743358, 0.821961896532945824, 0.697296697907720420,  1.095253089696273996, 1.694511991011381502, 1.996998717589755890, 2.273018402331320242, -0.191100605848157606, -0.191093614586216293, -0.294966549027094738, -0.607032364877862496, -0.849198688290699222, -0.974475621298580452, -1.375803201942292109, -1.290680675256194565, -0.786638229180541160, -0.398442469525911669,  0.402491553317749773 ))

  expect_equal(c(fit$state_vars),
               c(0.300276014967140625, 0.084325233765547583, 0.084325233765220997, 0.370338648573870444, 0.086354234348502298, 0.004087321675549315, 0.004087321675205677, 0.090747196638368249, 0.090904158884560313, 0.011852245581609001, 0.011852245581609017, 0.100108204393884487, 0.091783438022739153, 0.012412004658343193, 0.012412004658343200, 0.101236984080608577, 0.093031462387435843, 0.012651885439268214, 0.012651885439268169, 0.102598209056767170, 0.093203706005060130, 0.013226269212376587, 0.013226269212376560, 0.102161604952115792, 0.093745748578751009, 0.013816112543168771, 0.013816112543168769, 0.103658433694816504, 0.094125673709005653, 0.014203683868298962, 0.014203683868298959, 0.104945310257841068, 0.095009553824941015, 0.015590281291152964, 0.015590281291152952, 0.105335253029802317, 0.099815421709878463, 0.014885346087120999, 0.014885346087120999, 0.109742391438562542, 0.129845053223049289, 0.016736378070242752, 0.016736378070242766, 0.139987667950181360 ))

  expect_equal(c(fit$lag_one_cor),
               c(  8.444083390544754e-02, 3.255550270336604e-03, 3.264104801959131e-03, 8.812975361800844e-02, 3.213613609325543e-02, -3.703469736656911e-03,  -3.819864011710457e-03, 2.923214276873729e-02, 3.356879750491641e-02, -1.507987590652113e-03, -1.522676598773977e-03, 3.194701444439850e-02,   3.423140465138890e-02, -1.359467765867548e-03, -1.364599923285139e-03, 3.264045872400097e-02, 3.466405420013452e-02, -1.114686276046398e-03,  -1.189609798379931e-03, 3.287567959548075e-02, 3.482650174499644e-02, -8.233990692327227e-04, -7.940357498380981e-04, 3.312741801066402e-02,   3.493088077924721e-02, -5.118978495668373e-04, -4.739975991990813e-04, 3.383414873104589e-02, 3.431324496428288e-02, 4.463481355594556e-04,   3.481518685431802e-04, 3.368512731134585e-02, 3.002397989940960e-02, 1.797858797865795e-03, 1.816268713617769e-03, 3.046327421672968e-02,  1.978557419476790e-316, 6.952446316923950e-310, 6.951021612464943e-310, 6.952447816557452e-310 ))

  expect_equal(unname(c(fit$fixed_effects)),
               c(-3.6887661722538625, -0.3536299452646241 ))

  expect_equal(c(fit$order),
               c(1 ))

  expect_equal(c(fit$F_),
               c(1, 0, 0, 1 ))

  expect_equal(c(fit$method),
               c("UKF" ))

  expect_equal(c(fit$model),
               c("logit" ))

  expect_equal(c(fit$est_Q_0),
               c(FALSE ))

  expect_equal(c(fit$LR),
               c(1 ))

  # We just check that the calls succeeds
  expect_no_error(
    fit <- ddhazard(formula(survival::Surv(tstart, tstop, event) ~
                              -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + x2 + x3),
                    data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
                    control = list(method = "UKF")))


  # matplot(sims$betas, type = "l", ylim = range(fit$state_vecs, sims$betas))
  # matplot(fit$state_vecs, type = "l", col = 3:4, add = T, lty = 1)
  # abline(h = fit$fixed_effects, col = 1:2)
})

if(!interactive()) test_that("Mixtures versus previous fit for both types of models", expect_true(F))
if(!interactive()) test_that("Model with only one random", expect_true(F))
if(!interactive()) test_that("Manually check offsets in both parts of the algorithm", expect_true(F))
