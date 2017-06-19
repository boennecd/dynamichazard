context("Testing print function")

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
    a_0 = rep(0, 1), Q_0 = diag(1, 1),
    max_T = 20, order = 1)

  # print(paste0(capture.output(print(result, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(result, digits = 4),
              prints_text("Formula:\nsurvival::Surv(stop, event) ~ ddFixed(1) + group\n\nEstimated with EKF in 2 iterations of the EM algorithm\n\nEstimated time-varying effects and point-wise standard deviation:\n     group1    sd \nt0  -0.2848 0.7478\nt1  -0.3929 0.5696\nt2   0.1458 0.6336\nt3   0.8916 0.5785\nt4   0.6187 0.4406\nt5   1.7126 0.5019\nt6   1.5240 0.4089\nt7   0.6616 0.4473\nt8   0.8427 0.5509\nt9   0.6898 0.5493\nt10  0.5305 0.5821\nt11 -0.0174 0.6183\nt12 -0.1406 0.7227\nt13  0.3111 0.7925\nt14  1.5061 0.7855\nt15  1.0142 0.5459\nt16  0.3840 0.6165\nt17  0.1672 0.7149\nt18  0.4952 0.8014\nt19  0.6538 0.8563\nt20  1.5255 1.0509\n\nFixed effects are estimated in the E_step. The estimates are:\n(Intercept) \n     -2.925 ",
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
              prints_text("Bootstrap Statistics :\n                original     bias    bias (truncated)  std. error\n(Intercept):t0   -2.9663 -0.0215921        -0.0215921      0.2123\n(Intercept):t1   -2.9663 -0.0215923        -0.0215923      0.2123\n(Intercept):t2   -2.8843 -0.0124441        -0.0124441      0.2190\n(Intercept):t3   -2.7851 -0.0113302        -0.0113302      0.2128\n(Intercept):t4   -2.6954 -0.0124093        -0.0124093      0.2395\n(Intercept):t5   -2.6008 -0.0054419        -0.0054419      0.2639\n(Intercept):t6   -2.6131  0.0001238         0.0001238      0.2706\n(Intercept):t7   -2.7090  0.0029724         0.0029724      0.2590\n(Intercept):t8   -2.7866  0.0010536         0.0010536      0.2492\n(Intercept):t9   -2.8571  0.0043026         0.0043026      0.2388\n(Intercept):t10  -2.9262  0.0085269         0.0085269      0.2408\n(Intercept):t11  -3.0013  0.0053449         0.0053449      0.2307\n(Intercept):t12  -3.0564  0.0078671         0.0078671      0.2238\n(Intercept):t13  -3.0943  0.0060994         0.0060994      0.2166\n(Intercept):t14  -3.1005  0.0061054         0.0061054      0.2141\n(Intercept):t15  -3.1262 -0.0051858        -0.0051858      0.2074\n(Intercept):t16  -3.1572 -0.0078072        -0.0078072      0.2022\n(Intercept):t17  -3.1788 -0.0124789        -0.0124789      0.1958\n(Intercept):t18  -3.1764 -0.0169721        -0.0169721      0.1915\n(Intercept):t19  -3.1825 -0.0166572        -0.0166572      0.1906\n(Intercept):t20  -3.1685 -0.0159302        -0.0159302      0.1938\ngroup1:t0         0.4309  0.0426120         0.0426120      0.2422\ngroup1:t1         0.4309  0.0426118         0.0426118      0.2422\ngroup1:t2         0.4662  0.0489150         0.0489150      0.2454\ngroup1:t3         0.5160  0.0546436         0.0546436      0.2491\ngroup1:t4         0.5497  0.0553742         0.0553742      0.2414\ngroup1:t5         0.6037  0.0555932         0.0555932      0.2444\ngroup1:t6         0.6084  0.0567070         0.0567070      0.2358\ngroup1:t7         0.5729  0.0510402         0.0510402      0.2310\ngroup1:t8         0.5635  0.0458764         0.0458764      0.2276\ngroup1:t9         0.5480  0.0441578         0.0441578      0.2308\ngroup1:t10        0.5311  0.0439654         0.0439654      0.2278\ngroup1:t11        0.5090  0.0420390         0.0420390      0.2326\ngroup1:t12        0.5004  0.0419378         0.0419378      0.2389\ngroup1:t13        0.5041  0.0426926         0.0426926      0.2480\ngroup1:t14        0.5213  0.0451209         0.0451209      0.2595\ngroup1:t15        0.5172  0.0397781         0.0397781      0.2632\ngroup1:t16        0.5107  0.0373364         0.0373364      0.2683\ngroup1:t17        0.5124  0.0350343         0.0350343      0.2748\ngroup1:t18        0.5234  0.0332983         0.0332983      0.2834\ngroup1:t19        0.5309  0.0331911         0.0331911      0.2910\ngroup1:t20        0.5462  0.0337530         0.0337530      0.3010",
                          fixed = T))

  result <-  ddhazard(
    formula = survival::Surv(stop, event) ~ ddFixed(group),
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F, eps = 1e-1,
                   save_data = T, save_risk_set = T,
                   fixed_terms_method = "E_step"),
    a_0 = 0, Q_0 = 1e4, Q = 1e-2,
    max_T = 20, order = 1)

  set.seed(1992)
  suppressWarnings(boot_out <- ddhazard_boot(result, R = 20))

  # print(paste0(capture.output(print(boot_out, digits = 4)), collapse = "\n"), max.print = 1e8)
  expect_that(print(boot_out, digits = 4),
              prints_text("Bootstrap Statistics :\n                original   bias    bias (truncated)  std. error\n(Intercept):t0   -2.9719  0.02795           0.02795      0.2005\n(Intercept):t1   -2.9719  0.02795           0.02795      0.2005\n(Intercept):t2   -2.8892  0.04089           0.04089      0.2095\n(Intercept):t3   -2.7887  0.04762           0.04762      0.2260\n(Intercept):t4   -2.6945  0.05911           0.05911      0.2422\n(Intercept):t5   -2.6002  0.07923           0.07923      0.2670\n(Intercept):t6   -2.6097  0.07752           0.07752      0.2526\n(Intercept):t7   -2.7021  0.06988           0.06988      0.2366\n(Intercept):t8   -2.7803  0.06150           0.06150      0.2231\n(Intercept):t9   -2.8512  0.05553           0.05553      0.2116\n(Intercept):t10  -2.9209  0.04826           0.04826      0.2004\n(Intercept):t11  -2.9957  0.03659           0.03659      0.1967\n(Intercept):t12  -3.0515  0.02811           0.02811      0.1956\n(Intercept):t13  -3.0913  0.02760           0.02760      0.1950\n(Intercept):t14  -3.1005  0.02842           0.02842      0.1976\n(Intercept):t15  -3.1262  0.03090           0.03090      0.1926\n(Intercept):t16  -3.1567  0.03245           0.03245      0.1951\n(Intercept):t17  -3.1789  0.03214           0.03214      0.2008\n(Intercept):t18  -3.1779  0.03299           0.03299      0.2100\n(Intercept):t19  -3.1851  0.03352           0.03352      0.2168\n(Intercept):t20  -3.1730  0.03506           0.03506      0.2275\nddFixed(group)1   0.5247 -0.05077          -0.05077      0.2536",
                          fixed = T))
})
