context("Running tests for hatvalues")

#####
# Load data
if(interactive()){
  root <- gsub("(.+dynamichazard)(.*)", "\\1", getwd())
  diag_data_path <- paste0(root, "/vignettes/Diagnostics")
} else{
  diag_data_path <- "."
}

load(paste0(diag_data_path, "/Rossi.RData"))
load(paste0(diag_data_path, "/whas500.RData"))

#####
# Run tests

test_that("hatvalues works with dynamic effects only and gives previous results", {
  dd_fit <- ddhazard(
    Surv(start, stop, event) ~ fin + age + prio + employed.cumsum,
    data = Rossi, id = Rossi$id, by = 1, max_T = 52,
    Q_0 = diag(1000, 5), Q = diag(1e-2, 5),
    control = ddhazard_control(eps = 2e-3))

  # plot(dd_fit)

  hats <- hatvalues(dd_fit)
  hats <- hats[c(1, 3, 5)]
  # save_to_test(hats, "hats_dym_Rossi")

  expect_equal(hats, read_to_test("hats_dym_Rossi"), tolerance = 1.490116e-08)

  #####
  dd_fit <- ddhazard(
    Surv(lenfol, fstat) ~ gender + age + bmi + hr + cvd,
    data = whas500, by = 100, max_T = 2000,
    Q_0 = diag(10000, 6), Q = diag(.1, 6))

  # plot(dd_fit)

  hats <- hatvalues(dd_fit)
  hats <- hats[c(1, 3, 5)]
  # save_to_test(hats, "hats_dym_whas500")

  expect_equal(hats, read_to_test("hats_dym_whas500"), tolerance = 1.490116e-08)
})

test_that("hatvalues works with dynamic and fixed effects and gives previous results",{
  dd_fit <- ddhazard(
    Surv(start, stop, event) ~ ddFixed(fin) + age + prio + employed.cumsum,
    data = Rossi, id = Rossi$id, by = 1, max_T = 52,
    Q_0 = diag(1000, 4), Q = diag(.01, 4),
    control = ddhazard_control(fixed_terms_method = "E_step"))

  # dd_fit$fixed_effects
  # plot(dd_fit)

  hats <- hatvalues(dd_fit)
  hats <- hats[c(1, 3, 5)]
  # save_to_test(hats, "hats_dym_n_fixed_Rossi")

  expect_equal(hats, read_to_test("hats_dym_n_fixed_Rossi"), tolerance = 1.490116e-08)

  #####
  dd_fit <- ddhazard(
    Surv(lenfol, fstat) ~ ddFixed(gender) + ddFixed(age) + bmi + hr + cvd,
    data = whas500, by = 100, max_T = 1800,
    Q_0 = diag(10000, 4), Q = diag(.1, 4),
    control = ddhazard_control(fixed_terms_method =  'M_step', eps = 5e-3))

  # dd_fit$fixed_effects
  # plot(dd_fit)

  hats <- hatvalues(dd_fit)
  hats <- hats[c(1, 3, 5)]
  # save_to_test(hats, "hats_dym_n_fixed_whas500")

  expect_equal(hats, read_to_test("hats_dym_n_fixed_whas500"), tolerance = 1.490116e-08)
})
