context("Testing weights in fit")

test_that("Passing fewer weigths than rows in design mat throws error",{

  expect_error(
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1,
      control = ddhazard_control(
        est_Q_0 = F, n_max = 10^4, eps = 10^-4,
        save_data = F, save_risk_set = F),
      a_0 = rep(0, 2), Q_0 = diag(1, 2),
      max_T = 45,
      id = head_neck_cancer$id, order = 1,
      weights = rep(1, nrow(head_neck_cancer) - 1)
    )
  )
})

test_that("Using weights yields a message about lag_one_cov", {
  set.seed(11)
  ws <- sample.int(nrow(head_neck_cancer), 60)
  ws <- sapply(1:nrow(head_neck_cancer), function(x) sum(ws == x))
  ws[1] <- 2

  expect_message(
    ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 3,
      a_0 = rep(0, 2), Q_0 = diag(1, 2), Q = diag(1, 2),
      max_T = 27, order = 1,
      weights = ws),
    regexp = "^lag_one_cov will not be correct when some weights are not 1\n$")
})

test_that("Making large design mat and using weights yield the same with EKF",{
  set.seed(9191)
  tmp <- sample.int(nrow(head_neck_cancer), 100, replace = T)
  dum_design <- head_neck_cancer[tmp, ]

  ws <- sapply(1:nrow(head_neck_cancer), function(x) sum(tmp == x))

  meth <- "EKF"
  for(m in c("logit", exp_model_names)){
    f1 <- ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = dum_design,
      by = 5, model = m,
      a_0 = c(-2,0), Q_0 = diag(100, 2), Q = diag(1e-2, 2),
      control = ddhazard_control(method = meth, denom_term = 1e-3),
      max_T = 25, order = 1)

    suppressMessages(f2 <- ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 5, model = m,
      a_0 = c(-2, 0), Q_0 = diag(100, 2), Q = diag(1e-2, 2),
      control = ddhazard_control(method = meth, denom_term = 1e-3),
      max_T = 25, order = 1,
      weights = ws))

    info <- paste("m =", m)
    expect_equal(f1$state_vecs, f2$state_vecs, info = info, tolerance = 1e-5)
    expect_equal(f1$state_vars, f2$state_vars, info = info, tolerance = 1e-5)
  }
})


test_that("Making large design mat and using weights yield the same with UKF",{
  set.seed(9191)
  tmp <- sample.int(nrow(head_neck_cancer), 100, replace = T)
  dum_design <- head_neck_cancer[tmp, ]

  ws <- sapply(1:nrow(head_neck_cancer), function(x) sum(tmp == x))

  meth <- "UKF"
  for(m in c("logit", "exponential")){
    f1 <- ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = dum_design,
      by = 5, model = m,
      a_0 = c(-2,0), Q_0 = diag(1, 2), Q = diag(1e-2, 2),
      control = ddhazard_control(method = meth, denom_term = 1e-3),
      max_T = 25, order = 1)

    suppressMessages(f2 <- ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer,
      by = 5, model = m,
      a_0 = c(-2, 0), Q_0 = diag(1, 2), Q = diag(1e-2, 2),
      control = ddhazard_control(method = meth, denom_term = 1e-3),
      max_T = 25, order = 1,
      weights = ws))

    info <- paste("m =", m)
    expect_equal(f1$state_vecs, f2$state_vecs, info = info, tolerance = 1e-5)
    expect_equal(f1$state_vars, f2$state_vars, info = info, tolerance = 1e-5)
  }
})
