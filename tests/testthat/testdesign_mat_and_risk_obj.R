# library(survival); library(testthat); source("R/test_utils.R")
library(survival)

# Simulate data
set.seed(11111)
sims = as.data.frame(test_sim_func_logit(n_series = 10^5, n_vars = 10, beta_start = 1,
                                         intercept_start = - 16, sds = c(sqrt(.1), rep(1, 10)),
                                         x_range = 1, x_mean = .5)$res)


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


ids_to_row_indicies <- tapply(seq_len(nrow(sims)), sims$id, identity, simplify = F)
for(do_truncate in c(F, T)){
  # Truncate on second round
  if(do_truncate){
    sims$tstart = round(sims$tstart, digits = 2)
    sims$tstop = round(sims$tstop, digits = 2)

    sims = sims[!(sims$tstart == sims$tstop), ]

    # Have to get the design matrix again
    design = get_design_matrix(formula(Surv(tstart, tstop, event) ~ x1 + x2 +x3), sims)
  }

  for(by_ in c(1, 2.63)){
    # We start by testing without a max_T argument
    risk_set = get_risk_obj(Y = design$Y, by = by_, id = sims[, "id"])

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

        expect_true(all((design$Y[not_set, 1] > t_ - by_ | # after interval start
                            design$Y[not_set, 2] <= t_ - by_) | # stop befor interval start
                          (design$Y[not_set, 1] <= t_ - by_ & # is inside a bin
                             design$Y[not_set, 2] <= t_)))
        })
    }

    test_that("All risk objs is_event_in with no T_max",{
      expect_equal(sum(design$Y[, 3] == 1),
                   sum(risk_set$is_event_in > -1))
    })

    test_that("Testing risk sets default max_T and number of events", {
      expect_lte(max(design$Y[, 2][design$Y[, 3] == 1]), length(risk_set$risk_sets) * by_)
      expect_gte(max(design$Y[, 2][design$Y[, 3] == 1]), (length(risk_set$risk_sets) - 1) * by_)

      expect_equal(sims$event, design$Y[, 3], check.attributes = F, use.names = F)
      })

    inside_of_bin <- !seq_len(nrow(sims)) %in% is_in_a_bin
    test_that("All those inside of bins and that has the correct is_event_in", {
      expect_equal(risk_set$is_event_in[inside_of_bin], rep(-1, sum(inside_of_bin)))
    })

    inside_of_bin_and_event <- inside_of_bin & sims$event
    test_that("All those inside of bins with events have has the correct previous row as event", {
      ids <- sims$id[inside_of_bin_and_event]
      ids_to_check <- sample(ids, size = min(1e3, length(ids)), replace = F)
      for(id in ids_to_check){
        rows_ <- sims[ids_to_row_indicies[[id]], ]
        event_time <- max(rows_$tstop)
        bin_start <- min(risk_set$event_times[risk_set$event_times >= event_time]) - by_

        should_be_event <- which(rows_$tstart <= bin_start &  bin_start < rows_$tstop)

        if(risk_set$is_event_in[ids_to_row_indicies[[id]][should_be_event]] < 0)
          expect_true(FALSE)
      }
    })

    # Test with max_T
    max_T = 6
    max_T = (max_T - max_T %% by_) + by_ * (max_T %% by_ > 0)
    risk_set = get_risk_obj(Y = design$Y, by = by_, id = sims[, "id"], max_T = max_T)

    test_that("All risk objs is_event_in with no T_max",{
      expect_equal(sum(design$Y[, 3] == 1 & design$Y[, 2] <= max_T),
                   sum(risk_set$is_event_in > -1))
    })

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

        expect_true(all((design$Y[not_set, 1] > t_ - by_ | # after interval start
                           design$Y[not_set, 2] <= t_ - by_) | # stop befor interval start
                          (design$Y[not_set, 1] <= t_ - by_ & # is inside a bin
                             design$Y[not_set, 2] <= t_)))
      })
    }

    test_that("Testing the number of risk sets with max_T set",
              expect_equal(length(risk_set$risk_sets), max_T / by_))

    is_event = which(tapply(sims$tstop <= max_T & sims$event, sims$id, max) == 1) # is event within the max_T
    is_event_in_risk_set = which(tapply(risk_set$is_event_in > -1, sims$id, max) == 1) # is event within the max_T

    test_that("Those that have a failure with max_T also have a failure in the risk set", {
      expect_equal(is_event, is_event_in_risk_set)
    })

    inside_of_bin <- !seq_len(nrow(sims)) %in% is_in_a_bin & sims$tstop <= max_T
    test_that("All those inside of bins and that has the correct is_event_in", {
      expect_equal(risk_set$is_event_in[inside_of_bin], rep(-1, sum(inside_of_bin)))
    })

    inside_of_bin_and_event <- inside_of_bin & sims$event & sims$tstart < max_T
    test_that("All those inside of bins with events have has the correct previous row as event", {
      ids <- sims$id[inside_of_bin_and_event]
      ids_to_check <- sample(ids, size = min(1e3, length(ids)), replace = F)
      for(id in ids_to_check){
        rows_ <- sims[ids_to_row_indicies[[id]], ]
        event_time <- max(rows_$tstop)
        bin_start <- min(risk_set$event_times[risk_set$event_times >= event_time]) - by_

        should_be_event <- which(rows_$tstart <= bin_start &  bin_start < rows_$tstop)

        if(risk_set$is_event_in[ids_to_row_indicies[[id]][should_be_event]] < 0)
          expect_true(FALSE)
      }
    })

  }
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

test_that("is_for_discrete_model logic work",
          expect_true(FALSE))
