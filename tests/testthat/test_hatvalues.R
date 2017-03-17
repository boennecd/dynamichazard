# Had issues with win builder. Thus, these lines
test_name <- "hatvalues"
cat("\nRunning", test_name, "\n")

#####
# Load data
if(interactive()){
  diag_data_path <- paste0(stringr::str_extract(getwd(), ".+dynamichazard"), "/vignettes/Diagnostics")
} else{
  print(getwd()) # TODO: Delete
  diag_data_path <- "./vignettes/Diagnostics"
}

load(paste0(diag_data_path, "/Rossi.RData"))
load(paste0(diag_data_path, "/whas500.RData"))

#####
# Run tests

test_that("hatvalues works with dynamic effects only and gives previous results", {
  dd_fit <- ddhazard(
    Surv(start, stop, event) ~ fin + age + prio + employed.cumsum,
    data = Rossi, id = Rossi$id, by = 1, max_T = 52,
    Q_0 = diag(10000, 5), Q = diag(.1, 5))

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
    Q_0 = diag(10000, 4), Q = diag(.1, 4))

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
    control = list(fixed_terms_method =  'E_step'))

  # dd_fit$fixed_effects
  # plot(dd_fit)

  hats <- hatvalues(dd_fit)
  hats <- hats[c(1, 3, 5)]
  # save_to_test(hats, "hats_dym_n_fixed_whas500")

  expect_equal(hats, read_to_test("hats_dym_n_fixed_whas500"), tolerance = 1.490116e-08)
})

# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
