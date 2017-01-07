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
    a_0 = rep(0, 2), Q_0 = diag(1000, 2), Q = diag(1e-2, 2),
    max_T = 20, order = 1)

  set.seed(999)
  suppressWarnings(boot_out <- ddhazard_boot(result, R = 20))

  # print(paste0(capture.output(print(boot_out, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(boot_out, digits = 4),
              prints_text("Bootstrap Statistics :\n                original     bias    bias (truncated)  std. error\n(Intercept):t0   -2.9663 -0.0070150        -0.0070150      0.1907\n(Intercept):t1   -2.9663 -0.0070154        -0.0070154      0.1907\n(Intercept):t2   -2.8843 -0.0027237        -0.0027237      0.2040\n(Intercept):t3   -2.7851  0.0200975         0.0200975      0.2244\n(Intercept):t4   -2.6954  0.0249371         0.0249371      0.2464\n(Intercept):t5   -2.6008  0.0398451         0.0398451      0.2727\n(Intercept):t6   -2.6131  0.0310506         0.0310506      0.2737\n(Intercept):t7   -2.7090  0.0291893         0.0291893      0.2565\n(Intercept):t8   -2.7866  0.0256196         0.0256196      0.2339\n(Intercept):t9   -2.8571  0.0136321         0.0136321      0.2105\n(Intercept):t10  -2.9262  0.0113894         0.0113894      0.1908\n(Intercept):t11  -3.0013  0.0111960         0.0111960      0.1678\n(Intercept):t12  -3.0564  0.0129090         0.0129090      0.1476\n(Intercept):t13  -3.0943  0.0170695         0.0170695      0.1292\n(Intercept):t14  -3.1005  0.0228135         0.0228135      0.1163\n(Intercept):t15  -3.1262  0.0227104         0.0227104      0.1089\n(Intercept):t16  -3.1572  0.0164915         0.0164915      0.1073\n(Intercept):t17  -3.1788  0.0123075         0.0123075      0.1069\n(Intercept):t18  -3.1764  0.0083041         0.0083041      0.1133\n(Intercept):t19  -3.1825  0.0007634         0.0007634      0.1185\n(Intercept):t20  -3.1685 -0.0069382        -0.0069382      0.1294\ngroup1:t0         0.4309  0.0164983         0.0164983      0.2600\ngroup1:t1         0.4309  0.0164983         0.0164983      0.2600\ngroup1:t2         0.4662  0.0180187         0.0180187      0.2567\ngroup1:t3         0.5160  0.0245545         0.0245545      0.2540\ngroup1:t4         0.5497  0.0245076         0.0245076      0.2551\ngroup1:t5         0.6037  0.0257983         0.0257983      0.2564\ngroup1:t6         0.6084  0.0254617         0.0254617      0.2495\ngroup1:t7         0.5729  0.0258200         0.0258200      0.2506\ngroup1:t8         0.5635  0.0280522         0.0280522      0.2533\ngroup1:t9         0.5480  0.0247730         0.0247730      0.2548\ngroup1:t10        0.5311  0.0268654         0.0268654      0.2579\ngroup1:t11        0.5090  0.0287665         0.0287665      0.2634\ngroup1:t12        0.5004  0.0312556         0.0312556      0.2709\ngroup1:t13        0.5041  0.0341161         0.0341161      0.2805\ngroup1:t14        0.5213  0.0375155         0.0375155      0.2924\ngroup1:t15        0.5172  0.0367224         0.0367224      0.2999\ngroup1:t16        0.5107  0.0331945         0.0331945      0.3038\ngroup1:t17        0.5124  0.0294679         0.0294679      0.3089\ngroup1:t18        0.5234  0.0255338         0.0255338      0.3151\ngroup1:t19        0.5309  0.0209280         0.0209280      0.3211\ngroup1:t20        0.5462  0.0160353         0.0160353      0.3283",
                          fixed = T))

  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, n_max = 10^4, eps = 1e-1,
                   save_data = T, save_risk_set = T,
                   fixed_terms_method = "E_step"),
    a_0 = 0, Q_0 = 1e4, Q = 1e-2,
    max_T = 20, order = 1)

  set.seed(1992)
  suppressWarnings(boot_out <- ddhazard_boot(result, R = 20))

  # print(paste0(capture.output(print(boot_out, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(boot_out, digits = 4),
              prints_text("Bootstrap Statistics :\n                original    bias    bias (truncated)  std. error\n(Intercept):t0   -2.9727  0.013746          0.013746      0.2117\n(Intercept):t1   -2.9727  0.013746          0.013746      0.2117\n(Intercept):t2   -2.8895  0.025566          0.025566      0.2170\n(Intercept):t3   -2.7885  0.045818          0.045818      0.2254\n(Intercept):t4   -2.6938  0.059115          0.059115      0.2299\n(Intercept):t5   -2.5989  0.062517          0.062517      0.2243\n(Intercept):t6   -2.6085  0.064565          0.064565      0.2119\n(Intercept):t7   -2.7013  0.052130          0.052130      0.1976\n(Intercept):t8   -2.7799  0.038589          0.038589      0.1794\n(Intercept):t9   -2.8510  0.033308          0.033308      0.1678\n(Intercept):t10  -2.9210  0.031327          0.031327      0.1568\n(Intercept):t11  -2.9961  0.026648          0.026648      0.1462\n(Intercept):t12  -3.0521  0.022105          0.022105      0.1378\n(Intercept):t13  -3.0920  0.025620          0.025620      0.1375\n(Intercept):t14  -3.1012  0.031135          0.031135      0.1424\n(Intercept):t15  -3.1269  0.025482          0.025482      0.1429\n(Intercept):t16  -3.1575  0.015239          0.015239      0.1379\n(Intercept):t17  -3.1797  0.002558          0.002558      0.1321\n(Intercept):t18  -3.1787 -0.010517         -0.010517      0.1287\n(Intercept):t19  -3.1859 -0.012146         -0.012146      0.1262\n(Intercept):t20  -3.1737 -0.013756         -0.013756      0.1286\nddFixed(group)1   0.5246 -0.007686         -0.007686      0.2968",
                          fixed = T))
})






# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
