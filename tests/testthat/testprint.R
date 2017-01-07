if(interactive()){
  library(survival); library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")
}

# Had issues with win builder. Thus, these lines
test_name <- "print"
cat("\nRunning", test_name, "\n")

test_that("Print yields the expected results for object returned by ddhazard", {
  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    max_T = 20, order = 1)

  # print(paste0(capture.output(print(result, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(result, digits = 4),
              prints_text("Formula:\nsurvival::Surv(stop, event) ~ group\n\nEstimated with EKF in 2 iterations of the EM algorithm\n\nEstimated time-varying effects and point-wise standard deviation:\n    (Intercept)    sd   group1    sd \nt0       -1.802 0.6526 0.04333 0.7070\nt1       -2.337 0.2934 0.08234 0.3937\nt2       -2.803 0.3770 0.04737 0.4966\nt3       -2.785 0.4345 0.47529 0.5713\nt4       -2.434 0.4579 0.03807 0.5567\nt5       -1.785 0.4263 0.79952 0.5653\nt6       -1.973 0.3739 0.66402 0.4720\nt7       -2.450 0.4074 0.19807 0.5189\nt8       -2.878 0.4578 0.70811 0.6143\nt9       -3.016 0.5186 0.73993 0.6517\nt10      -3.094 0.5527 0.72452 0.6901\nt11      -3.389 0.5746 0.40006 0.7262\nt12      -3.599 0.6131 0.44561 0.8123\nt13      -3.733 0.6429 0.98993 0.8732\nt14      -3.116 0.6851 2.13339 0.8809\nt15      -3.279 0.6218 1.61627 0.7177\nt16      -3.576 0.6465 1.13830 0.7640\nt17      -3.894 0.6828 1.04528 0.8481\nt18      -3.689 0.7515 1.36305 0.9472\nt19      -3.950 0.7812 1.55307 1.0131\nt20      -3.613 0.9426 2.34269 1.2122",
                          fixed = T))

  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ ddFixed(1) + group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = F, save_risk_set = F,
                   fixed_terms_method = "E_step"),
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    max_T = 20, order = 1)

  # print(paste0(capture.output(print(result, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(result, digits = 4),
              prints_text("Formula:\nsurvival::Surv(stop, event) ~ ddFixed(1) + group\n\nEstimated with EKF in 3 iterations of the EM algorithm\n\nEstimated time-varying effects and point-wise standard deviation:\n      group2    sd  group1    sd \nt0  -0.12734 0.8745 0.1273 0.8745\nt1  -0.13878 0.9878 0.1450 1.0187\nt2   0.44852 1.0907 0.7447 1.1319\nt3   0.61420 1.1119 1.4339 1.1391\nt4   1.36766 1.1345 1.1935 1.0988\nt5   1.51739 1.1048 2.2492 1.1257\nt6   1.42554 1.1021 2.0740 1.0910\nt7   1.17621 1.1118 1.2379 1.1050\nt8   0.61468 1.1244 1.3992 1.1475\nt9   0.46213 1.1554 1.2493 1.1482\nt10  0.37433 1.1716 1.0877 1.1632\nt11  0.29088 1.1817 0.5633 1.1810\nt12  0.15197 1.1893 0.4445 1.2336\nt13 -0.12769 1.1994 0.8584 1.2714\nt14 -0.18141 1.2264 1.9502 1.2675\nt15  0.03043 1.2439 1.5042 1.1525\nt16  0.03152 1.2341 0.9260 1.1841\nt17 -0.16720 1.2300 0.7390 1.2339\nt18 -0.15531 1.2582 1.0564 1.2826\nt19 -0.40867 1.2899 1.2180 1.3189\nt20 -0.54259 1.3846 2.0253 1.4447\n\nFixed effects are estimated in the E_step. The estimates are:\n(Intercept) \n     -3.489 ",
                          fixed = T))

  suppressWarnings(result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ ddFixed(1) + ddFixed(as.numeric(group)),
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = F, save_risk_set = F,
                   fixed_terms_method = "E_step"),
    max_T = 20, order = 1))

  # print(paste0(capture.output(print(result, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(result, digits = 4),
              prints_text("Formula:\nsurvival::Surv(stop, event) ~ ddFixed(1) + ddFixed(as.numeric(group))\n\nEstimated with EKF in 3 iterations of the EM algorithm\n\nFixed effects are estimated in the E_step. The estimates are:\n      (Intercept) as.numeric(group) \n          -3.4628            0.5524 ",
                          fixed = T))
})

test_that("print.ddhazard_boot gives the expected output", {
  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = T, save_risk_set = T),
    a_0 = rep(0, 2), Q_0 = diag(1, 2),
    max_T = 20, order = 1)

  set.seed(999)
  suppressWarnings(boot_out <- ddhazard_boot(result, R = 20))

  # print(paste0(capture.output(print(boot_out, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(boot_out, digits = 4),
              prints_text("Bootstrap Statistics :\n                original    bias    bias (truncated)  std. error\n(Intercept):t0  -3.10368 -0.178792         -0.178792      0.3313\n(Intercept):t1  -3.30198 -0.211916         -0.211916      0.3361\n(Intercept):t2  -3.09440 -0.152364         -0.152364      0.4571\n(Intercept):t3  -2.78856  0.175317          0.175317      0.2815\n(Intercept):t4  -2.41954  0.229656          0.229656      0.6212\n(Intercept):t5  -1.91191  0.250943          0.250943      0.5484\n(Intercept):t6  -2.04226  0.061962          0.061962      0.3717\n(Intercept):t7  -2.45615  0.126424          0.126424      0.3816\n(Intercept):t8  -2.81828  0.028341          0.028341      0.3164\n(Intercept):t9  -2.97114 -0.097480         -0.097480      0.3226\n(Intercept):t10 -3.08504 -0.080145         -0.080145      0.2923\n(Intercept):t11 -3.32848 -0.040542         -0.040542      0.2734\n(Intercept):t12 -3.50014 -0.089217         -0.089217      0.3359\n(Intercept):t13 -3.61147 -0.240672         -0.240672      0.4456\n(Intercept):t14 -3.30824  0.104003          0.104003      0.4750\n(Intercept):t15 -3.37985  0.146607          0.146607      0.6705\n(Intercept):t16 -3.58288  0.109803          0.109803      0.6411\n(Intercept):t17 -3.79062  0.045142          0.045142      0.4580\n(Intercept):t18 -3.68967  0.082536          0.082536      0.5423\n(Intercept):t19 -3.83830 -0.049035         -0.049035      0.4220\n(Intercept):t20 -3.65394 -0.260416         -0.260416      0.5031\ngroup1:t0        0.05441  0.041273          0.041273      0.3099\ngroup1:t1        0.07776 -0.024491         -0.024491      0.6269\ngroup1:t2        0.18784  0.127172          0.127172      0.4834\ngroup1:t3        0.56202  0.419658          0.419658      0.5030\ngroup1:t4        0.12697 -0.077222         -0.077222      0.6754\ngroup1:t5        0.75512 -0.242025         -0.242025      1.0386\ngroup1:t6        0.67084 -0.054930         -0.054930      0.5000\ngroup1:t7        0.24083 -0.201424         -0.201424      0.4898\ngroup1:t8        0.64454  0.306500          0.306500      0.6970\ngroup1:t9        0.68299  0.125766          0.125766      0.5090\ngroup1:t10       0.66915  0.200462          0.200462      0.3982\ngroup1:t11       0.39442  0.008934          0.008934      0.4137\ngroup1:t12       0.43573 -0.143888         -0.143888      0.6595\ngroup1:t13       0.86854 -0.034538         -0.034538      0.5086\ngroup1:t14       1.72627 15.845453         15.845453     30.2872\ngroup1:t15       1.39674 15.705899         15.705899     30.2990\ngroup1:t16       1.03752 15.751224         15.751224     30.2192\ngroup1:t17       1.00113 15.156919         15.156919     29.0263\ngroup1:t18       1.28822 15.544379         15.544379     29.6715\ngroup1:t19       1.48486 14.585662         14.585662     27.9770\ngroup1:t20       2.10541 13.716398         13.716398     26.7153",
                          fixed = T))

  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = T, save_risk_set = T),
    a_0 = 0, Q_0 = 1,
    max_T = 20, order = 1)

  set.seed(1992)
  suppressWarnings(boot_out <- ddhazard_boot(result, R = 20))

  # print(paste0(capture.output(print(boot_out, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(boot_out, digits = 4),
              prints_text("Bootstrap Statistics :\n                original   bias    bias (truncated)  std. error\n(Intercept):t0   -3.5322 -0.07839          -0.07839      0.2365\n(Intercept):t1   -3.5312 -0.10276          -0.10276      0.2570\n(Intercept):t2   -3.0545  0.39813           0.39813      1.4237\n(Intercept):t3   -2.6628  0.26084           0.26084      0.6413\n(Intercept):t4   -2.3848  0.20276           0.20276      0.4452\n(Intercept):t5   -1.8065  0.36386           0.36386      0.9431\n(Intercept):t6   -1.9070  0.26248           0.26248      0.3884\n(Intercept):t7   -2.3841  0.16237           0.16237      0.2893\n(Intercept):t8   -2.6307  0.04388           0.04388      0.2250\n(Intercept):t9   -2.7641  0.04771           0.04771      0.2564\n(Intercept):t10  -2.8985  0.25765           0.25765      0.4024\n(Intercept):t11  -3.1567  0.15087           0.15087      0.2806\n(Intercept):t12  -3.2909  0.20223           0.20223      0.6286\n(Intercept):t13  -3.3296  0.21722           0.21722      0.3869\n(Intercept):t14  -3.0592  0.49790           0.49790      0.7323\n(Intercept):t15  -3.1041  0.37879           0.37879      0.5564\n(Intercept):t16  -3.2504  0.24275           0.24275      0.4500\n(Intercept):t17  -3.3538  0.03572           0.03572      0.3068\n(Intercept):t18  -3.2146 -0.05767          -0.05767      0.3241\n(Intercept):t19  -3.2489 -0.08269          -0.08269      0.3791\n(Intercept):t20  -3.0546  0.33471           0.33471      1.8794\nddFixed(group)1   0.2546 -0.37961          -0.37961      0.8565",
                          fixed = T))
})






# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
