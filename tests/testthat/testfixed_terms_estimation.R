library(biglm)

if(interactive()){
  library(survival); library(dynamichazard); library(testthat)
  source("C:/Users/boennecd/Dropbox/skole_backup/phd/dynamichazard/R/test_utils.R")
}

set.seed(548237)
sims <- test_sim_func_logit(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
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


test_that("Get previous results with logit model with some fixed terms", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + x1 + x2 + ddFixed(x3))

  suppressMessages(
    res1 <- ddhazard(form, data = sims$res, model = "logit", by = 1, id = sims$res$id, max_T = 10,
                     control = list(save_risk_set = F, save_data = F)))

  # matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
  # matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
  # abline(h = res1$fixed_effects, col = c(1,4))

  # get_expect_equal(res1, file = "tmp.txt")

  expect_equal(c(res1$state_vecs),
               c( 0.778541338485961631, 0.778879145699752140, 0.599450494552086099, 0.106689722725345881, 0.279157387889645225,  0.677148028899600463, 1.122303854092461695, 2.672609292005364789, 5.673221383772331983, 5.409865322374935559,  8.265469678693211364, -1.362247025595830552, -1.362749500199533559, -0.426376628773305844, 0.419392535713738412,  1.364398199564436753, 1.299657532244780977, 2.697040200314659852, 1.765037913448112006, 1.243299110123040618,  1.368462134119203810, -0.009284232971245698 ))

  expect_equal(c(res1$state_vars),
               c( 2.320121369669100986, -0.555196891151226946, -0.555196891151226946, 1.256137952305543593, 0.330138597842687675,  0.001129604695150169, 0.001129604695150173, 0.281146454737927098, 0.306509478720435491, -0.013702605999991013, -0.013702605999991010, 0.240325897541568911, 0.359534503146291540, -0.040524258123700869, -0.040524258123700869,  0.279134015350714249, 0.383192816080901499, -0.047037710157792650, -0.047037710157792643, 0.277587786732228003,  0.347135636509643408, -0.051921043419431191, -0.051921043419431191, 0.268504509620995646, 0.363939442182868889, -0.065381560065589650, -0.065381560065589650, 0.283024598415983930, 0.336013856780990772, -0.061144404823326046, -0.061144404823326046, 0.237147901951189122, 0.307146046159263442, -0.077878412973142921, -0.077878412973142921,  0.254147852817652231, 0.173595524943991997, -0.035009470890098815, -0.035009470890098815, 0.198886375795202713,  0.257870763004976522, -0.051706503748950826, -0.051706503748950826, 0.280431632079972792 ))

  expect_equal(c(res1$lag_one_cov),
               c( 0.2604225003571093433, 0.0173671275315504595, 0.0203608365149278783, 0.2522938778486077394,  0.0434757066594110564, 0.0207240781908055519, 0.0208036137708710565, 0.0646364914239786759,  0.0443939671086375151, 0.0155069451747530309, 0.0132201880091016889, 0.0619474168108883616,  0.0534298946905332819, 0.0106270538368620007, 0.0114845151497209863, 0.0691940395181235568,  0.0508830012709687862, 0.0092442956671195415, 0.0071300658249742205, 0.0660895133706539267,  0.0474212877544047395, 0.0050726960763465842, 0.0037661609582322190, 0.0665212710933604234,  0.0455260921494123685, 0.0013175369494681731, 0.0022837186600320543, 0.0584530234363593235,  0.0376887479836751291, -0.0002247849275741078, -0.0036761685770627214, 0.0520772615606443570,  0.0195653482731466544, -0.0010920190402388999, -0.0008629781404866463, 0.0446801906315306419,  0.0166742945713191162, 0.0013137733802308866, 0.0016138358552887635, 0.0504832801149985300 ))

  expect_equal(unname(c(res1$fixed_effects)),
               c(-3.843074021376860, -2.891431151165306  ))

  expect_equal(c(res1$n_iter),
               c(13 ))

  expect_equal(c(res1$Q),
               c( 2.7445708074529653, -0.8349012935426495, -0.8349012935426495, 1.2100343381078322 ))

  expect_equal(c(res1$Q_0),
               c(10, 0, 0, 10 ))

  expect_equal(c(res1$n_risk),
               c(1000, 972, 948, 933, 912, 887, 846, 810, 745, 693 ))

  expect_equal(c(res1$times),
               c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ))

  expect_equal(c(res1$risk_set),
               c(NULL ))

  expect_equal(c(res1$data),
               c(NULL ))

  expect_equal(c(res1$order),
               c(1 ))

  expect_equal(c(res1$F_),
               c(1, 0, 0, 1 ))

  expect_equal(c(res1$method),
               c("EKF" ))

  expect_equal(c(res1$model),
               c("logit" ))

  expect_equal(c(res1$est_Q_0),
               c(FALSE ))

  expect_equal(c(res1$LR),
               c(1 ))
})







set.seed(312237)
sims <- test_sim_func_exp(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
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

test_that("Get previous results with exponential model with some fixed terms", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + x1 + x2 + x3)

  suppressMessages(
    res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
                    control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3,
                                   save_risk_set = F, LR = .1, save_data = F, n_max = 1e3),
                    Q_0 = diag(rep(10, 3))))

  # matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
  # matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
  # abline(h = res1$fixed_effects, col = 1)
  # get_expect_equal(res1, file = "tmp.txt", eps = 1e-3)

  expect_equal(c(res1$state_vecs),
               c(-0.46536788127469964, -0.46534249186861132, 0.31311070660717083, 0.02649174129956280, -0.17168995901195203, -0.47579045473084852, -0.42184080213975428, -0.19505356116579509, 0.18645304725442680, 0.27516210574912919, -0.03710365743442967, -2.38851748160244837, -2.38852086381647277, -2.43893653388935139, -2.44497366155019158, -2.36900002774638585, -2.34697958348962121, -2.24518116682354263, -2.13042285761170458, -2.09666891358525076, -2.10848354031333507, -2.13556958391085772, -1.78254086986964277, -1.78253073513642279, -3.27980312620167647, -2.49670986234400827, -1.84236503399558815, -1.08859500258539188, -0.92134103892936503, -1.11396586895128902, -1.56822455411768802, -1.42399376205381034, -0.71053243765412311 )
               , tolerance = 0.001)

  expect_equal(c(res1$state_vars),
               c( 0.3880744834173412983, 0.0093948301702217376, -0.3869128803496821867, 0.0093948301702217445,  0.1532439775895504397, -0.0062558191689357032, -0.3869128803496822422, -0.0062558191689356911,  0.9916790339869070436, 0.1669551685334740854, 0.0037713073307852757, -0.0716408311757816968,  0.0037713073307852826, 0.1068044523904205789, -0.0168332785783618655, -0.0716408311757817107, -0.0168332785783618516, 0.2695256723595723658, 0.1346608840605048596, 0.0007487024846483913, -0.0670650852146195686, 0.0007487024846483947, 0.0857587463778033410, -0.0125575910329147125, -0.0670650852146195686, -0.0125575910329147195, 0.2319788703467001167, 0.1060676006335572741,  0.0050051398127107126, -0.0399716179875695340, 0.0050051398127107126, 0.0745643133676427150, -0.0189164608992849030, -0.0399716179875695271, -0.0189164608992848995, 0.1677502514094211250,  0.1116118430743910789, 0.0053304729396792981, -0.0564790249349131149, 0.0053304729396792989,  0.0717133537398907256, -0.0178516195169398373, -0.0564790249349131218, -0.0178516195169398373,  0.1933085291560894536, 0.1145222088943739502, 0.0047846892870662110, -0.0609019706478392750,  0.0047846892870662110, 0.0726333913705557027, -0.0131785450274266268, -0.0609019706478392750, -0.0131785450274266285, 0.2084220512444648343, 0.1285805351937768515, 0.0031206752647873056, -0.0778037996066952287, 0.0031206752647873065, 0.0754645372081586274, -0.0085846398420904423, -0.0778037996066952287, -0.0085846398420904423, 0.2409752356733285139, 0.1378912277895883443,  0.0032372679533255126, -0.0827077236383596537, 0.0032372679533255126, 0.0813053673170728020, -0.0074277132067602624, -0.0827077236383596537, -0.0074277132067602615, 0.2586899548280685024,  0.1470218695724367364, 0.0039120499999022815, -0.0784607459977116478, 0.0039120499999022815,  0.0901179155709007640, -0.0075359504996076821, -0.0784607459977116617, -0.0075359504996076812,  0.2607873417555857665, 0.1604358830413227799, 0.0071470797298289427, -0.0697926942907377973,  0.0071470797298289410, 0.1044401466976822845, -0.0121275637461476377, -0.0697926942907377834, -0.0121275637461476359, 0.2618834569783081712, 0.1996982887584216326, 0.0121629890335292999, -0.0811719039657386177, 0.0121629890335292999, 0.1300894426154544370, -0.0164459485150391591, -0.0811719039657386177, -0.0164459485150391591, 0.3280341563442024655 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$fixed_effects)),
               c(-3.803338792043399  )
               , tolerance = 0.001)

  expect_equal(c(res1$n_iter),
               c(21 )
               , tolerance = 0.001)

  expect_equal(c(res1$Q),
               c( 0.249097019901144801, 0.006356048073213742, -0.372582991975295152, 0.006356048073213742, 0.045169703702657216,  0.011158141930496983, -0.372582991975295152, 0.011158141930496983, 0.844304111946077462 )
               , tolerance = 0.001)

  expect_equal(c(res1$Q_0),
               c(10, 0, 0, 0, 10, 0, 0, 0, 10 )
               , tolerance = 0.001)

  expect_equal(c(res1$n_risk),
               c(2005, 1927, 1850, 1746, 1749, 1649, 1623, 1581, 1543, 1541 )
               , tolerance = 0.001)

  expect_equal(c(res1$times),
               c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 )
               , tolerance = 0.001)

  expect_equal(c(res1$risk_set),
               c(NULL )
               , tolerance = 0.001)

  expect_equal(c(res1$data),
               c(NULL )
               , tolerance = 0.001)

  expect_equal(c(res1$order),
               c(1 )
               , tolerance = 0.001)

  expect_equal(c(res1$F_),
               c(1, 0, 0, 0, 1, 0, 0, 0, 1 )
               , tolerance = 0.001)

  expect_equal(c(res1$method),
               c("EKF" )
               , tolerance = 0.001)

  expect_equal(c(res1$model),
               c("exponential" )
               , tolerance = 0.001)

  expect_equal(c(res1$est_Q_0),
               c(FALSE )
               , tolerance = 0.001)

  expect_equal(c(res1$LR),
               c(0.04444444444444445 )
               , tolerance = 0.001)
})

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
