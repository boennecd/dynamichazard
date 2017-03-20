if(interactive()){
  library(survival); library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("./R/test_utils.R")

  get_frac_n_weights <- with(environment(ddhazard), get_frac_n_weights)
}

# Had issues with win builder. Thus, these lines
test_name <- "boot"
cat("\nRunning", test_name, "\n")

test_that("Throws error when risk_set or data is not saved",{
  for(tmp in list(c(T, F),
                  c(F, T),
                  c(F, F))){
    fit <- ddhazard(
      formula = survival::Surv(stop, event) ~ group,
      data = head_neck_cancer, max_T = 40,
      by = 1, a_0 = c(0, 0),
      control = list(save_data = tmp[1], save_risk_set = tmp[2]))

    expect_error(ddhazard_boot(fit, 2),
                 regexp = "^Cannot bootstrap estimates when ddhazard has been")
  }
})

test_that("boot yields previously computed values with pbc", {
  suppressMessages(
    fit <- ddhazard(Surv(tstart, tstop, death == 2) ~ age + edema +
                      log(albumin) + log(protime) + log(bili), pbc2,
                    id = pbc2$id, by = 100, max_T = 3600,
                    Q_0 = diag(rep(10000, 6)), Q = diag(rep(0.001, 6)),
                    control = list(save_risk_set = T, save_data = T, eps = .1)))

  #####
  set.seed(993)
  suppressWarnings(tmp <- ddhazard_boot(fit, do_sample_weights = F, do_stratify_with_event = F, R = 99))

  expect_no_error(plot(fit, ddhazard_boot = tmp))

  tmp <- tmp[c("t0", "t")]
  tmp$t <- tmp$t[1:5, 1:100]
  # save_to_test(tmp, "boot1")

  expect_equal(tmp, read_to_test("boot1"))

  #####
  set.seed(994)
  suppressWarnings(tmp <- ddhazard_boot(fit, do_sample_weights = T, do_stratify_with_event = F, R = 99))

  expect_no_error(plot(fit, ddhazard_boot = tmp))

  tmp <- tmp[c("t0", "t")]
  tmp$t <- tmp$t[1:5, 1:100]
  # save_to_test(tmp, "boot2")

  expect_equal(tmp, read_to_test("boot2"))

  #####
  set.seed(995)
  suppressWarnings(tmp <- ddhazard_boot(fit, do_sample_weights = F, do_stratify_with_event = T, R = 99))

  expect_no_error(plot(fit, ddhazard_boot = tmp))

  tmp <- tmp[c("t0", "t")]
  tmp$t <- tmp$t[1:5, 1:100]
  # save_to_test(tmp, "boot3")

  expect_equal(tmp, read_to_test("boot3"))

  #####
  set.seed(999)
  suppressWarnings(tmp <- ddhazard_boot(fit, do_sample_weights = T, do_stratify_with_event = T, R = 99))

  expect_no_error(plot(fit, ddhazard_boot = tmp))

  tmp <- tmp[c("t0", "t")]
  tmp$t <- tmp$t[1:5, 1:100]
  # save_to_test(tmp, "boot4")

  expect_equal(tmp, read_to_test("boot4"))
})

######
# Test get_frac_n_weights function

test_that("frac_n_weights gives error if sample is too small", {
  expect_error(get_frac_n_weights(49, .99),
               "^Sample of 49 is too small to give 0.99 confidence bounds$")

  expect_error(get_frac_n_weights(49, .01),
               "^Sample of 49 is too small to give 0.01 confidence bounds$")
})

tmp <- get_frac_n_weights(999, .9519)
test_that("Weights are as expected", {
  R <- 991 - 1 # Large-ish prime number less one

  for(i in 1:99){
    a <- i / 100
    frac_n_weights <- get_frac_n_weights(R, a)
    lbl <- paste("i =", i)

    expect_lte(frac_n_weights$k, a * (R + 1), label = lbl)
    expect_gte(frac_n_weights$k + 1, a * (R + 1), label = lbl)
    expect_equal(frac_n_weights$w_k + frac_n_weights$w_k_p_1, 1, label = lbl)
    expect_true(all(c(frac_n_weights$w_k, frac_n_weights$w_k_p_1) > 0),
                label = lbl)
  }
})

test_that("Boot works with posterior_approx",
          expect_true(FALSE))

# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
