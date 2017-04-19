# Had issues with win builder. Thus, these lines
test_name <- "test_GMA"
cat("\nRunning", test_name, "\n")
options(ddhazard_use_speedglm = F)

if(interactive()){
  library(survival); library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("R/test_utils.R")
}

# frm <-  ~ age + edema + log(albumin) + log(bili) + log(protime)
# args <- list(formula = update(frm, Surv(tstart, tstop, death == 2) ~ .),
#              data = pbc2,
#              id = pbc2$id, by = 100, max_T = 3600,
#              control = list(method = "EKF"),
#              Q_0 = diag(rep(1e6, 6)), Q = diag(rep(1e-4, 6)))
#
# fit <- do.call(ddhazard, args)
# max(apply(fit$state_vars, 3, kappa))
#
# pca <- prcomp(frm, pbc2)
# pca$sdev
#
# pca_dat <- predict(pca, model.frame(frm, pbc2))
# pca_dat <- cbind(pbc2[, c("tstart", "tstop", "death", "id")], pca_dat)
#
# args$formula <- Surv(tstart, tstop, death == 2) ~ PC1 + PC2 + PC3
# args$data <- pca_dat
# args$id <- pca_dat$id
# args$Q_0 <- diag(rep(1e6, 4))
# args$Q <- diag(rep(1e-4, 4))
# args$control$method <- "GMA"
# args$control$debug <- TRUE
#
# fit <- do.call(ddhazard, args)
# max(apply(fit$state_vars, 3, kappa))

# tmp <- static_glm(Surv(tstart, tstop, death == 2) ~ age + edema +
#                     log(albumin) + log(bili) + log(protime), pbc2,
#                   id = pbc2$id, by = 100, max_T = 3600)
# kappa(vcov(tmp))
# qr(vcov(tmp))
# car::vif(tmp)

test_that("Gives same results w/ 1. order logit model", {
  tmp <- file("tmp.txt")
  sink(tmp)
  fit <- do.call(ddhazard, args)
  sink()
  close(tmp)
})

test_that("GMA gives the same w/ all exponential model inputs", {
  expect_true(FALSE)
})

test_that("GMA works w/ second order random walk", {
  expect_true(FALSE)
})

test_that("GMA works w/ fixed effects in E-step", {
  expect_true(FALSE)
})

test_that("GMA works w/ fixed effects in M-step", {
  expect_true(FALSE)
})

test_that("Changing hyper parameters w/ GMA changes the result", {
  expect_true(FALSE)
})

test_that("GMA makes one mesasge when global scoring did not succed within given number of iterations", {
  expect_true(FALSE)
})


# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
