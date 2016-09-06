# library(survival); library(testthat); source("R/test_utils.R")
library(survival)

# Simulate data
set.seed(11111)
sims = as.data.frame(test_sim_func_logit(n_series = 10^5, n_vars = 10, beta_start = 1,
                                         intercept_start = -10, sds = c(sqrt(.1), rep(1, 10)))$res)



#######
# Test design mat
design = get_design_matrix(formula(Surv(tstart, tstop, event) ~ x1 + x2 +x3), sims)

test_that("Testing get design marix", {
  expect_equal(design$Y[, 1], sims$tstart, check.attributes = F, use.names = F)
  expect_equal(design$Y[, 2], sims$tstop, check.attributes = F, use.names = F)
  expect_equal(design$Y[, 3], sims$event, check.attributes = F, use.names = F)

  expect_equal(colnames(design$X), c("(Intercept)", "x1", "x2", "x3"))
  expect_equal(design$X[, 2:4], as.matrix(sims[, c("x1", "x2", "x3")]), check.attributes = F)
  })


############
# test rcpp get risk set function
by_ = 1

for(do_truncate in c(F, T)){
  # Truncate on second round
  if(do_truncate){
    sims$tstart = round(sims$tstart, digits = 2)
    sims$tstop = round(sims$tstop, digits = 2)

    sims = sims[!(sims$tstart == sims$tstop), ]

    # Have to get the design matrix again
    design = get_design_matrix(formula(Surv(tstart, tstop, event) ~ x1 + x2 +x3), sims)
  }

  # We start by testing without a max_T argument
  risk_set = get_risk_obj(Y = design$Y, by = by_, id = sims[, "id"])

  t_ = 0
  for(set in risk_set$risk_sets){
    t_ = t_ + by_

    set = unlist(set)
    not_set = setdiff(seq_len(nrow(sims)), set)

    test_that("Testing risk sets with no max_T", {
      expect_true(all(design$Y[set, 2] > t_ - by_))
      expect_true(all(design$Y[set, 1] < t_))

      expect_true(all(design$Y[not_set, 1] > t_ - by_ | # after interval start
                        design$Y[not_set, 2] <= t_ - by_)) # stop befor interval start
      })
  }

  test_that("Testing risk sets default max_T and number of events", {
    expect_lt(max(design$Y[, 2][design$Y[, 3] == 1]), length(risk_set$risk_sets) * by_)
    expect_gte(max(design$Y[, 2][design$Y[, 3] == 1]), (length(risk_set$risk_sets) - 1) * by_)

    expect_equal(sims$event, design$Y[, 3], check.attributes = F, use.names = F)
    })

  # Test with max_T
  max_T = 3
  risk_set = get_risk_obj(Y = design$Y, by = by_, id = sims[, "id"], max_T = max_T)

  t_ = 0
  is_in_a_bin = c()
  for(set in risk_set$risk_sets){
    t_ = t_ + by_

    set = unlist(set)
    not_set = setdiff(seq_len(nrow(sims)), set)

    is_in_a_bin = union(set, is_in_a_bin)

    test_that("Testing risk sets with no max_T", {
      expect_true(all(design$Y[set, 2] > t_ - by_))
      expect_true(all(design$Y[set, 1] < t_))

      expect_true(all(design$Y[not_set, 1] > t_ - by_ | # start after interval start
                        design$Y[not_set, 2] <= t_ - by_)) # stop befor interval start
    })
  }

  test_that("Testing the number of risk sets with max_T set",
            expect_equal(length(risk_set$risk_sets), 3))

  is_event = which(tapply(sims$tstop < max_T & sims$event, sims$id, max) == 1) # is event within the max_T

  test_that("Testing that new flags are valid when max_T is set",{
    is_new_flag = which(risk_set$new_events_flags & !sims$event)
    for(f in is_new_flag){
      id = sims$id[f]
      id_rows = which(sims$id == id)

      if(!risk_set$new_events_flags[max(intersect(is_in_a_bin, id_rows))])
        expect_true(FALSE)
    }

    expect_equal(sum(!risk_set$new_events_flags & sims$event), 0)
  })
}


#######
# Further test for design mat

# Test without Y
arg_list <- list(formula = formula(Surv(tstart, tstop, event) ~ x1 + x2 + x3), data = sims, response = T)
design <- do.call(get_design_matrix, arg_list)

arg_list$response <- F
design_no_Y <- do.call(get_design_matrix, arg_list)

test_that("get_design_matrix with versus without response", {
  expect_equal(design$X, design_no_Y$X)
  expect_false(design$formula == design_no_Y$formula)
  expect_null(design_no_Y$Y)
})

arg_list$data <- sims[, !colnames(sims) %in% c("tstart", "tstop", "event")]

tryCatch({
  test_that("calling get_design_matrix without response and also no response in data",{
    design_no_Y_also_in_data <- do.call(get_design_matrix, arg_list)

    expect_equal(design_no_Y_also_in_data$X, design_no_Y$X)
    expect_equal(design_no_Y_also_in_data$formula, design_no_Y$formula)
    expect_null(design_no_Y_also_in_data$Y)
  })
}, error = function(...){
  test_that("calling get_design_matrix without response and also no response in data",
            expect_true(FALSE))
})

# Test that removal of intercepts works as in lm
arg_list <- list(formula = formula(Surv(tstart, tstop, event) ~ x1 + x2 + x3), data = sims, response = T)
design <- do.call(get_design_matrix, arg_list)

arg_list$formula <- formula(Surv(tstart, tstop, event) ~ -1 + x1 + x2 + x3)
design_no_intercept <- do.call(get_design_matrix, arg_list)

test_that("Design mat with vs. without intercept", {
  expect_equal(design$Y, design_no_intercept$Y)
  expect_equal(design$X[, 2:4], design_no_intercept$X[, 1:3])
  expect_equal(ncol(design_no_intercept$X[, 1:3]), 3)
})
