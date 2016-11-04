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

  # get_expect_equal(fit, file = "get_expect_equal.txt")

  expect_equal(names(fit$fixed_effects), c("rep(1, length(x1))", "x1"))

  expect_equal(c(fit$state_vecs),
               c( 0.11865797440743088, 0.11866022563897899, 0.21224697794330311,  0.51328205431702356, 0.82960127689169250, 0.98838046698727589,  1.04155102454529258, 1.46285200407928295, 1.93214959313360679,  2.15778831672109073, 2.31254709793364466, -0.32751677529358308, -0.32752390699433587, -0.48306055239490253, -0.73697550252305999, -0.91372277657383194, -0.98947233003915747, -1.11352113687604737, -0.88074571214538366, -0.42600145174710302, -0.05881452883697125,  0.44550170865199346 ))

  expect_equal(c(fit$state_vars),
               c(0.56499571961569117, 0.17453009657509888, 0.17453009657519883, 0.62236703803853466, 0.27826343876836690, 0.07785789639592139, 0.07785789639602814, 0.30415158848966029, 0.28442534267867131, 0.08284967044928834, 0.08284967044928834, 0.31103515429504913, 0.28555897886892667, 0.08308004951709749, 0.08308004951709745, 0.31202939143161207, 0.28678788179201020, 0.08330918132208770, 0.08330918132208767, 0.31319942684891711, 0.28747893076708608, 0.08389542864891690, 0.08389542864891686, 0.31373968875956215, 0.28966947846461855, 0.08448137410382769, 0.08448137410382769, 0.31666468328881053, 0.29476234848579741, 0.08554533122543628, 0.08554533122543628, 0.32261543753717831, 0.30863047872557048, 0.08925482375908786, 0.08925482375908786, 0.33685482024372559, 0.34801319702082339, 0.09866433830980131, 0.09866433830980129, 0.37875741242326666, 0.45888611643030702, 0.13221361567744483, 0.13221361567744483, 0.49950417162687261 ))

  expect_equal(c(fit$lag_one_cor),
               c( 2.688970867191910e-01, 7.222062491087741e-02, 7.222608288655824e-02,  2.929326181277456e-01, 1.647879366118103e-01, 4.073830637393284e-02,  4.070767600968590e-02, 1.781152341063062e-01, 1.678415658592686e-01,  4.299172583878002e-02, 4.298508860619533e-02, 1.814154135129742e-01,  1.685278765348665e-01, 4.321225258639694e-02, 4.320827304882893e-02,  1.820268625344066e-01, 1.683925457095060e-01, 4.359742990735148e-02,  4.358626617220988e-02, 1.818216439777128e-01, 1.669182538671503e-01,  4.399302937560495e-02, 4.400882985250514e-02, 1.805980008293575e-01,  1.624051850878922e-01, 4.390768944001947e-02, 4.392333701485661e-02,  1.765749065809095e-01, 1.489309115116473e-01, 4.225257315108288e-02,  4.224102465634497e-02, 1.626007811461578e-01, 1.103830134450838e-01,  3.341782952451835e-02, 3.341716400444748e-02, 1.210649701220078e-01, 2.022171163176555e-316, 6.951897076961135e-310, 0.000000000000000e+00, 6.951899172369823e-310 ))

  expect_equal(unname(fit$fixed_effects),
               c(-3.6986970275194344, -0.3541344883005885 ))

  expect_equal(c(fit$n_iter),
               c(23 ))

  expect_equal(c(fit$Q),
               c(0.3157256648647212, 0.1152333387591277, 0.1152333387591277, 0.3522474006503304 ))

  # matplot(sims$betas, type = "l", ylim = range(fit$state_vecs, sims$betas))
  # matplot(fit$state_vecs, type = "l", col = 3:4, add = T, lty = 1)
  # abline(h = fit$fixed_effects, col = 1:2)

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
