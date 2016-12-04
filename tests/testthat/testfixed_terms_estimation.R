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
               c( 0.774026566564674945, 0.774366869252967938, 0.597980492507363293, 0.106360547493080215, 0.280336566432623113,  0.677663713014745817, 1.122299115291607086, 2.668452487869186651, 5.675589064223387936, 5.409825400516144356,  8.272148432063712420, -1.359356994089573112, -1.359860977427649331, -0.423828858843989575, 0.423365461990504621,  1.363075866039706430, 1.300909706816193001, 2.699600628983378936, 1.775221622693042622, 1.244895659894105222,  1.369194537603316331, -0.007574960829139332 ))

  expect_equal(c(res1$state_vars),
               c( 2.325808035942558583, -0.556314056954815506, -0.556314056954815506, 1.255898438469721157, 0.331467086299748470,  0.001005572779382430, 0.001005572779382419, 0.281998085260910136, 0.307702397495087765, -0.013885161830733577, -0.013885161830733577, 0.240908549432121449, 0.360999562577125843, -0.040817702124379650, -0.040817702124379650,  0.279807492691421800, 0.384739052381862734, -0.047340938493736641, -0.047340938493736634, 0.278202665133104432,  0.348455734508277493, -0.052217258615905773, -0.052217258615905773, 0.269127440045229982, 0.365343541307478259, -0.065745079402884640, -0.065745079402884654, 0.283657568068083565, 0.337140385954186939, -0.061460603953837192, -0.061460603953837192, 0.237559057533206663, 0.308209002310561864, -0.078413992015304870, -0.078413992015304870,  0.254542167316882395, 0.173832631284923433, -0.035130204020244506, -0.035130204020244506, 0.199088463663849330,  0.258395246081691399, -0.051913734398399414, -0.051913734398399414, 0.280900237025662025 ))

  expect_equal(c(res1$lag_one_cov),
               c( 0.2612928845783098475, 0.0173521922047875365, 0.0203668760661175723, 0.2530781259746648915, 0.0436811386780537603,  0.0208204151963864020, 0.0209006777578030038, 0.0650344986182617990, 0.0446023826017280245, 0.0155677431080261802,  0.0132619162308012951, 0.0622937521630072197, 0.0536884205274128015, 0.0106533150380225955, 0.0115179155097433726,  0.0695645726614575188, 0.0511161851446121987, 0.0092729057235389056, 0.0071377665571012014, 0.0664404109660966524,  0.0476292280545493885, 0.0050740036778336604, 0.0037554002496294862, 0.0668759525610890387, 0.0457050307790657570,  0.0012951530840262099, 0.0022625256432185045, 0.0587329271320696580, 0.0378155737559588909, -0.0002748790134440952, -0.0037620029414713934, 0.0522883142379131916, 0.0196103481820070341, -0.0011555033481836722, -0.0009052822251955746,  0.0448712918291203772, 0.0192595565122363060, 0.0003162259821144417, 0.0002281551077286639, 0.0519400029975638014 ))

  expect_equal(unname(c(res1$fixed_effects)),
               c(-3.843624231834462, -2.892056257728873  ))

  expect_equal(c(res1$n_iter),
               c(13 ))

  expect_equal(c(res1$Q),
               c( 2.753200630645960, -0.837098919178177, -0.837098919178177, 1.209077219027612 ))

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
               c( 0.05229798341008220, 0.05230081544449205, 0.06549627388696208, 0.07062882235727586, 0.06832559408518291,  0.07067667700685010, 0.08654655606312064, 0.09497545919567900, 0.09916793978431987, 0.09569154385177298,  0.08972938430407998, -3.77435992644363916, -3.77437180650135673, -3.78093563778847708, -3.78073098273714114, -3.76569813775232021, -3.75905629174024103, -3.73647543953614303, -3.71250131964646579, -3.70290182084688535, -3.70018112058471704, -3.70453144661074285, -3.00730868970687437, -3.00728470728469421, -3.02240054635580391, -3.00939169321823163, -2.97246731788115603, -2.92549473136960003, -2.89109170499457324, -2.86244202767785172, -2.83989965446351178, -2.81566503807353596, -2.79555796417696323 )
               , tolerance = 0.001)

  expect_equal(c(res1$state_vars),
               c( 0.1455000301664846063, 0.0013872101315510340, 0.0020285388502197959, 0.0013872101315510342, 0.1666738232535305286, -0.0960001725150314428, 0.0020285388502197942, -0.0960001725150316093, 0.2036566790382270398, 0.1208716539572260440,  0.0006619820062629626, 0.0017286885288622673, 0.0006619820062629626, 0.1389475940174205171, -0.0841047501882874471,  0.0017286885288622655, -0.0841047501882876136, 0.1673010393982573563, 0.1046785816974427907, 0.0003812176032880997,  0.0014144617794020955, 0.0003812176032880997, 0.1220532513086110826, -0.0758847105283917689, 0.0014144617794020964, -0.0758847105283918244, 0.1445386456084940541, 0.0944494213245384295, 0.0006010504346533273, 0.0012647092057449351,  0.0006010504346533265, 0.1118032946310972275, -0.0703699850054126663, 0.0012647092057449343, -0.0703699850054126802,  0.1303734971241726492, 0.0884863747900483943, 0.0008274990749969425, 0.0006705192250375772, 0.0008274990749969420,  0.1065181516085665497, -0.0672963204979765139, 0.0006705192250375776, -0.0672963204979765278, 0.1216633040602481258,  0.0863523292259017056, 0.0013408533297532815, 0.0003946752669257594, 0.0013408533297532815, 0.1051464734056778011, -0.0661491247401871652, 0.0003946752669257594, -0.0661491247401871652, 0.1184309288107977948, 0.0881230721637346615,  0.0016555917377680139, -0.0001577108994320903, 0.0016555917377680139, 0.1072080674607605588, -0.0672418496059406234, -0.0001577108994320905, -0.0672418496059406096, 0.1204716235352210202, 0.0931624998772399848, 0.0021417197557732004, -0.0003756519150862247, 0.0021417197557732004, 0.1130159214476939150, -0.0704137405419832524, -0.0003756519150862249, -0.0704137405419832385, 0.1270809607596524449, 0.1014061949107363247, 0.0024690756187012728, -0.0001155717103049341,  0.0024690756187012723, 0.1222399990768931466, -0.0753359900175574609, -0.0001155717103049341, -0.0753359900175574609,  0.1380588373003872049, 0.1134711705269920434, 0.0027879176038267468, 0.0003532153083166324, 0.0027879176038267468,  0.1361199689964802118, -0.0826924930442165668, 0.0003532153083166324, -0.0826924930442165668, 0.1556086724585228898,  0.1308894049675259075, 0.0029158643776528662, 0.0008041390876047570, 0.0029158643776528662, 0.1564536184718270384, -0.0932152310818315788, 0.0008041390876047570, -0.0932152310818315788, 0.1820491100853802169 )
               , tolerance = 0.001)

  expect_equal(c(res1$lag_one_cov),
               c(NULL )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$fixed_effects)),
               c(-4.291795701752389)
               , tolerance = 0.001)

  expect_equal(c(res1$n_iter),
               c(26 )
               , tolerance = 0.001)

  expect_equal(c(res1$Q),
               c( 0.0241115828599779494, 0.0007326253822521967, 0.0003407405379612195, 0.0007326253822521967, 0.0275570253088919248, -0.0124093484212931017, 0.0003407405379612195, -0.0124093484212931017, 0.0364984021890618973 )
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
               c(0.1 )
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
                  control = list(method = "UKF", fixed_parems_start = rep(0, 2),
                                 save_data = F, save_risk_set = F))


  # matplot(sims$betas, type = "l", lty = 1)
  # matplot(fit$state_vecs, type = "l", lty = 2, col = 3:4, add = T)
  # abline(h = fit$fixed_effects, col = 1:2)
  # get_expect_equal(fit, file = "tmp.txt")

  expect_equal(c(fit$state_vecs),
               c( 0.007448614508951475, 0.007449296611761227, -0.008343730452707976, 0.275988484153517621,  0.650607726330696745, 0.822035933445071310, 0.697722147350732125, 1.095577594153456458,  1.694614723669126022, 1.997196091824424258, 2.273328971932273124, -0.191047483635335413, -0.191040473807200128, -0.294991170064621000, -0.607010251510704535, -0.849166346484256929, -0.974543274949183802, -1.375828554842977836, -1.290695781316861135, -0.786763830718261947, -0.398495118139918614, 0.402284357656369829 ))

  expect_equal(c(fit$state_vars),
               c(0.300327982791047887, 0.084321286611075991, 0.084321286610892027, 0.370422456036598291, 0.086413500713397873, 0.004091096370294816, 0.004091096370101248, 0.090812701657028386, 0.090962198518812432, 0.011860458938372759, 0.011860458938179193, 0.100177581550918687, 0.091841808529020075, 0.012421553102273920, 0.012421553102080367, 0.101307687629529436, 0.093089770330813670, 0.012661760300479722, 0.012661760300286193, 0.102669172464966624, 0.093261482584170219, 0.013236541693320631, 0.013236541693127072, 0.102231958844678200, 0.093802985105780418, 0.013827488666673295, 0.013827488666479707, 0.103728539365780734, 0.094181347374331287, 0.014216289299019061, 0.014216289298825431, 0.105014014537682526, 0.095065378573909315, 0.015603872547911311, 0.015603872547717636, 0.105402541819223178, 0.099878009176012666, 0.014899556060138618, 0.014899556059944962, 0.109814860482982002, 0.129930782266901579, 0.016756540609913773, 0.016756540609720080, 0.140086188194331274 ))

  expect_equal(c(fit$lag_one_cov),
               c(NULL ))

  expect_equal(unname(fit$fixed_effects),
               c(-3.6887808556085830, -0.3536289109191706  ))

  expect_equal(c(fit$n_iter),
               c(12 ))

  expect_equal(c(fit$Q),
               c(0.2236200548757658, 0.0865644965834774, 0.0865644965834774, 0.2952452769482359 ))

  expect_equal(c(fit$Q_0),
               c(10, 0, 0, 10 ))

  expect_equal(c(fit$n_risk),
               c(10000, 9832, 9666, 9480, 9272, 9039, 8709, 8433, 8151, 7880 ))

  expect_equal(c(fit$times),
               c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ))

  expect_equal(c(fit$risk_set),
               c(NULL ))

  expect_equal(c(fit$data),
               c(NULL ))

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
