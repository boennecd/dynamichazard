library(biglm)

if(interactive()){
  library(survival); library(dynamichazard); library(testthat)
  source("C:/Users/boennecd/Dropbox/skole_backup/phd/dynamichazard/R/test_utils.R")
}

# Had issues with win builder. Thus, these lines
test_name <- "testfixed_terms_est_in_M_step"
cat("\nRunning", test_name, "\n")

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
                     control = list(eps_fixed_parems = 1e-3, fixed_effect_chunk_size = 1e3, max_it_fixed_params = 10)))

  tmp_design <- get_survival_case_weigths_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10)

  suppressWarnings(res2 <- bigglm(update(form, Y ~ .), data = tmp_design$X, family = binomial(), chunksize = 1e3))

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

  # get_expect_equal(res1, file = "tmp.txt", )

  expect_equal(c(res1$state_vecs),
               c( 1.33784167880023874, 1.33823452523718789, 1.35674083984908522, 0.24908802502559804, 0.49373599633551574,  0.80959321901894421, 1.93598617766876102, 2.22229643721847703, 5.78018898380291102, 5.32528177435792127,  8.29142611104554561, -0.94192897273469389, -0.94230230954427607, -0.42036367259209273, 0.31349922346842579,  1.21432275031608072, 1.20075200900525592, 2.38778033747193019, 2.28131864555201025, 0.06869842095379919,  1.36612737407374585, 0.76361111448619878 ))

  expect_equal(c(res1$state_vars),
               c( 2.5527695148530451519, -0.7162607232060363982, -0.7162607232060363982, 1.4697629370302198737, 0.3332722051483297254, -0.0002652750557076977, -0.0002652750557076977, 0.3049573913874989439, 0.3152337745513269507, -0.0172135290748056807, -0.0172135290748056807, 0.2638165216981178673, 0.3188671472130887863, -0.0326990784700045625, -0.0326990784700045556,  0.2725164877454505641, 0.3962889482328874613, -0.0703219773737278897, -0.0703219773737278897, 0.2964686997596947537,  0.3585297941396810018, -0.0627184542281121749, -0.0627184542281121749, 0.2874841771479983987, 0.3641340973890704324, -0.0611189361430429165, -0.0611189361430429096, 0.2942602028363880962, 0.3262257454021991188, -0.0795034693070020837, -0.0795034693070020837, 0.2596420956520920642, 0.3265389620231562984, -0.1002383460774469714, -0.1002383460774469714,  0.2743856275875475870, 0.2074903361612033537, -0.0338072612570374173, -0.0338072612570374173, 0.2233501160179429124,  0.2817995470735897290, -0.0614878578714110821, -0.0614878578714110821, 0.3154370280842789809 ))

  expect_equal(c(res1$lag_one_cov),
               c( 0.2547217836312979933, 0.0228375593958701044, 0.0249470320873037762, 0.2673198466264695705, 0.0412571139410628507,  0.0230703017320416301, 0.0234076834200742240, 0.0661886273480801479, 0.0372050173604432238, 0.0173218022357521859,  0.0153838533643861179, 0.0575388222696163465, 0.0435071132467827609, 0.0111804879078115024, 0.0109956336101321769,  0.0611806164518808732, 0.0474104499145256500, 0.0076852842281753229, 0.0066243400676007777, 0.0628844414959601611,  0.0438012880517647182, 0.0077148731495179021, 0.0078475962703599585, 0.0629263182718022424, 0.0387272148809793687,  0.0043394516254513359, 0.0019901758384756035, 0.0558602543635311199, 0.0334633900495338935, -0.0026803852757651099, -0.0050173815986381725, 0.0504869562098130692, 0.0216036239953910740, -0.0006726407844395997, 0.0017652858577329811,  0.0450870472087231758, 0.0230909365391175243, 0.0044745858767311899, 0.0020357878895393565, 0.0559590808289664732 ))

  expect_equal(unname(c(res1$fixed_effects)),
               c(-3.898278628172867, -3.151152030185528  ))

  expect_equal(c(res1$n_iter),
               c(15 ))

  expect_equal(c(res1$Q),
               c( 3.196659336700227, -1.141232138074276, -1.141232138074276, 1.515680568932240 ))

  expect_equal(c(res1$Q_0),
               c(10, 0, 0, 10 ))

  expect_equal(c(res1$n_risk),
               c(1000, 978, 957, 937, 913, 887, 844, 811, 751, 707 ))

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







set.seed(4682146)
sims <- test_sim_func_exp(n_series = 1e3, n_vars = 3, t_0 = 0, t_max = 10,
                          x_range = 1, x_mean = 0, re_draw = T, beta_start = 0,
                          intercept_start = -4, sds = c(.1, rep(1, 3)))
# sum(sims$res$event)

test_that("Only fixed effects yields same results as bigglm with exponential model", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + ddFixed(x1) + ddFixed(x2) + ddFixed(x3))

  suppressWarnings(res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
                                    control = list(eps_fixed_parems = 1e-4, fixed_effect_chunk_size = 1e3)))

  tmp_design <- get_survival_case_weigths_and_data(form, data = sims$res, by = 1, id = sims$res$id,
                                                   use_weights = F, max_T = 10, is_for_discrete_model = F)

  suppressWarnings(res2 <- bigglm(update(form, Y ~ . + offset(log(pmin(tstop, t) - pmax(tstart, t - 1)))),
                                  data = tmp_design$X, family = poisson(), chunksize = 1e3,
                                  tolerance = 1e-4))

  expect_equal(unname(coef(res2)), unname(c(res1$fixed_effects))
               , tolerance = 1e-05)
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
  arg_list_tmp$control$max_it_fixed_params <- 5
  suppressWarnings(res2 <- do.call(ddhazard, arg_list_tmp))
  expect_true(!all(res2$fixed_effects == res1$fixed_effects))
})

test_that("Get previous results with exponential model with some fixed terms", {
  form <- formula(survival::Surv(tstart, tstop, event) ~
                    -1 + ddFixed(rep(1, length(x1))) + x1 + x2 + x3)

  suppressMessages(
    res1 <- ddhazard(form, data = sims$res, model = "exponential", by = 1, id = sims$res$id, max_T = 10,
                    control = list(eps_fixed_parems = 1e-12, fixed_effect_chunk_size = 1e3,
                                   save_risk_set = F, save_data = F, n_max = 1e2,
                                   ridge_eps = .005),
                    Q_0 = diag(rep(10, 3)), Q = diag(.1, 3)))

  # matplot(sims$betas, type = "l", ylim = range(sims$betas, res1$state_vecs))
  # matplot(res1$state_vecs, add = T, col = 2:4, type = "l", lty = 1)
  # abline(h = res1$fixed_effects, col = 1)
  # get_expect_equal(res1, file = "tmp.txt", eps = 1e-3)

  expect_equal(unname(c(res1$state_vecs)),
               c(-1.2409255967381210, -1.2410119952430347, -0.8823363340358323, -1.8930694299170967, -1.3869475808530671, -1.5596828615517266, -1.3882636605713645, -2.8070572464401478, -4.4473231688135249, -3.8593189000088710, -3.7087893506239045, 1.4048635329562620, 1.4049813283630117, 0.9756419739592664, 2.3258771881941818,  1.5439836643687870, 1.6460768330695259, 1.3411527569106574, 3.1411996734552625, 5.2506354127311266,  4.5082492716541491, 4.4073869458674020, -1.5881448741926907, -1.5881809970940601, -1.5402411480601206, -1.7658041779958307, -1.5037516970722784, -1.4354883623068331, -1.3050650441198479, -1.5879550623486005, -1.9441409895787478, -1.9163631318139993, -1.9907085157313082 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$state_vars)),
               c( 1.100499998013223e+00, -1.099519896161443e+00, 1.073101898268469e-01, -1.099519896161443e+00,  1.646649238366633e+00, -2.420971406099759e-01, 1.073101898268469e-01, -2.420971406099759e-01,  3.278705345434609e-01, 3.804146580438273e-01, -2.180962552145562e-01, -2.545817391020100e-02, -2.180962552145562e-01, 4.696862283648261e-01, -3.428125732927978e-02, -2.545817391020097e-02, -3.428125732927981e-02, 2.328993139619872e-01, 3.261513994664113e-01, -1.896321961328797e-01, -1.279555770585577e-02, -1.896321961328797e-01, 4.193501488909193e-01, -3.255965039094157e-02, -1.279555770585575e-02, -3.255965039094158e-02, 1.965583110343035e-01, 3.845439321634586e-01, -2.846683033748543e-01, 9.928797450641106e-03, -2.846683033748543e-01, 5.218401550897300e-01, -4.967884289849327e-02, 9.928797450641120e-03, -4.967884289849328e-02, 1.755342363759886e-01,  2.505005296332327e-01, -1.249248911567435e-01, -2.428276858382328e-02, -1.249248911567435e-01,  3.004075746133845e-01, 3.915583028465300e-03, -2.428276858382327e-02, 3.915583028465307e-03,  1.456814425977130e-01, 3.168954722937026e-01, -2.196755485864360e-01, -3.380339618450356e-04, -2.196755485864360e-01, 4.227508229375952e-01, -2.745683593849374e-02, -3.380339618450252e-04, -2.745683593849374e-02, 1.453954673600227e-01, 3.069571487100030e-01, -2.117589811035098e-01,  3.517820265049952e-03, -2.117589811035098e-01, 4.147333044467880e-01, -3.285689455671033e-02,  3.517820265049983e-03, -3.285689455671031e-02, 1.415471926999000e-01, 3.522648341123036e-01, -2.737058459253242e-01, 1.891503936165169e-02, -2.737058459253242e-01, 5.008240087642825e-01, -5.230678734315113e-02, 1.891503936165168e-02, -5.230678734315112e-02, 1.398908114072103e-01,  2.410205715713883e-01, -1.327432871515942e-01, 8.977083710520056e-05, -1.327432871515942e-01,  3.260741710988039e-01, -2.407335598320027e-02, 8.977083710519362e-05, -2.407335598320028e-02,  1.261036844188326e-01, 1.133649092057718e-01, 3.369638959296151e-02, -1.467284390039982e-02,  3.369638959296148e-02, 1.209859620498711e-01, 4.717762792367312e-04, -1.467284390039983e-02,  4.717762792367321e-04, 1.140040843607593e-01, 1.700421939373131e-01, -3.086722604408289e-03, -2.480685386094261e-02, -3.086722604408289e-03, 1.821851550693598e-01, -2.011591212245160e-03, -2.480685386094261e-02, -2.011591212245160e-03, 1.439414079777564e-01 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$lag_one_cov)),
               c(NULL )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$fixed_effects)),
               c(-4.451438078513003 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$n_iter)),
               c(25 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$Q)),
               c( 1.0286096666960212, -1.2775628011667084, 0.1877614804143118, -1.2775628011667084, 1.6891270728429437, -0.2815121477803785, 0.1877614804143118, -0.2815121477803785, 0.1080165222675747 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$Q_0)),
               c(10, 0, 0, 0, 10, 0, 0, 0, 10 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$n_risk)),
               c(1956, 1931, 1882, 1856, 1808, 1758, 1771, 1706, 1722, 1608 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$times)),
               c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$risk_set)),
               c(NULL )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$data)),
               c(NULL )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$order)),
               c(1 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$F_)),
               c(1, 0, 0, 0, 1, 0, 0, 0, 1 )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$method)),
               c("EKF" )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$model)),
               c("exponential" )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$est_Q_0)),
               c(FALSE )
               , tolerance = 0.001)

  expect_equal(unname(c(res1$LR)),
               c(0.1975308641975309 )
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
               c(-0.022767706972769645, -0.022783587830587308, 0.006570684894700973, 0.615621587505349455, 1.329374104221104247,  1.196460676446490456, 0.731106826542401556, 1.128580235783753993, 1.861441882714902185, 1.829569953498984214,  1.772984046725965657, -0.088918260521527112, -0.088928522893182316, -0.554180842809645302, -0.693216624212211086, -1.107683977690919530, -0.997906043675983789, -1.644962214520343746, -0.849763339181839661, -0.749442170019782972, -0.203034923456535593, 0.955010801900676665 ))

  expect_equal(c(fit$state_vars),
               c(0.3024564385141008671, 0.0241394694996310138, 0.0241394695018546483, 0.4456378978995783058, 0.0435828059935903997, 0.0010172021281515226, 0.0010172021305311970, 0.0440386897362189736, 0.0397734542955975601, 0.0010414029575227602, 0.0010414029599023140, 0.0431095868104318033, 0.0400679505215062598, 0.0008470591314661429, 0.0008470591338456257, 0.0433086116432563964, 0.0399706078048868221, 0.0010589294903258187, 0.0010589294927053499, 0.0438036284035986587, 0.0372048711134163462, 0.0029852266982740641, 0.0029852267006536463, 0.0410303287719140652, 0.0391400819526037044, 0.0024000951529275097, 0.0024000951553070585, 0.0435961335146979184, 0.0407046780111072204, 0.0013631033208520377, 0.0013631033232316592, 0.0427529361148437995, 0.0427859053687857996, 0.0017722917205026379, 0.0017722917228822931, 0.0467355970340707985, 0.0403674790166447875, 0.0028672358206676164, 0.0028672358230472524, 0.0460823084065290420, 0.0473176083320425156, 0.0011078262379435530, 0.0011078262403230767, 0.0538300387648574796 ))

  expect_equal(c(fit$lag_one_cov),
               c(NULL ))

  expect_equal(unname(c(fit$fixed_effects)),
               c(-3.7285423369293724, -0.4070337753961206  ))

  expect_equal(c(fit$n_iter),
               c(8 ))

  expect_equal(c(fit$Q),
               c(0.26984202181816397, 0.02523298628396323, 0.02523298628396323, 0.42410431455626368 ))

  expect_equal(c(fit$Q_0),
               c(10, 0, 0, 10 ))

  expect_equal(c(fit$n_risk),
               c(10000, 9846, 9668, 9501, 9258, 9010, 8695, 8452, 8197, 7926 ))

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






# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
