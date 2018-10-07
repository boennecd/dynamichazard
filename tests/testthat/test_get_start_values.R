context("Testing 'get_start_values'")

test_that("'get_start_values' gives the same as 'static_glm' when it should with differnt families", {
  exs <- list(
    binomial = list(dt = quote(logit_sim_200), fam = "logit"),
    expo     = list(dt = quote(exp_sim_200)  , fam = "exponential"))

  for(e in exs){
    scl <- bquote(static_glm(
      formula = Surv(tstart, tstop, event) ~ x1 + x2, data = .(e$dt)$res,
      by = 1L, max_T = 10, family = .(e$fam), id = .(e$dt)$res$id))
    fit <- eval(scl)

    tmp <- eval(bquote(get_design_matrix_and_risk_obj(
      data = .(e$dt)$res, by = 1L, max_T = 10L,
      is_for_discrete_model = switch (e$fam,
        logit = TRUE, exponential = FALSE, stop()),
      id = .(e$dt)$res$id, fixed = Surv(tstart, tstop, event) ~ x1,
      random = ~ x2 - 1)))

    cl <- bquote(get_start_values(
      data = .(e$dt)$res, X = t(tmp$X_Y$X), max_T = 10, risk_set =
        tmp$risk_set, fixed_terms = t(tmp$X_Y$fixed_terms), order = 1,
      n_threads = 1L,
      fixed = Surv(tstart, tstop, event) ~ x1, random = ~ x2 - 1, type = "RW",
      model = .(e$fam)))
    r <- eval(cl)
    expect_equal(fit$coefficients, c(r$fixed_parems_start, r$a_0),
                 check.attributes = FALSE, info = e$fam)

    cl$a_0 <- 99
    r <- eval(cl)
    expect_equal(fit$coefficients[1:2], r$fixed_parems_start,
                 check.attributes = FALSE, info = e$fam)
    expect_equal(r$a_0, 99, info = e$fam)

    cl$a_0 <- NULL
    cl$fixed_parems_start <- c(99, 99)
    r <- eval(cl)
    expect_equal(c(99, 99), r$fixed_parems_start,
                 check.attributes = FALSE, info = e$fam)
    expect_equal(r$a_0, fit$coefficients[3], info = e$fam)

    scl$formula <- quote(Surv(tstart, tstop, event) ~ x1)
    cl$fixed_parems_start <- NULL
    cl$type <- "VAR"
    r <- eval(cl)
    fit <- eval(scl)
    expect_equal(fit$coefficients, r$fixed_parems_start,
                 check.attributes = FALSE, info = e$fam)
    expect_equal(r$a_0, 0, info = e$fam)

    cl$fixed_parems_start <- c(99, 99)
    r <- eval(cl)
    expect_equal(c(99, 99), r$fixed_parems_start,
                 check.attributes = FALSE, info = e$fam)

    cl$fixed_parems_start <- NULL
    cl$a_0 <- 99
    r <- eval(cl)
    expect_equal(fit$coefficients, r$fixed_parems_start,
                 check.attributes = FALSE, info = e$fam)
    expect_equal(r$a_0, 99, info = e$fam)
  }
})
