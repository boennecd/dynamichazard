if(interactive()){
  library(survival); library(testthat); source("R/test_utils.R"); library(dynamichazard)
}

get_design_matrix <- function(...) environment(ddhazard)$get_design_matrix(...)

# Simulate data
set.seed(11111)
sims = as.data.frame(test_sim_func_logit(n_series = 2e3, n_vars = 10, beta_start = 1,
                                         intercept_start = - 16, sds = c(sqrt(.1), rep(1, 10)),
                                         x_range = 1, x_mean = .5)$res)


#######
# Test design mat
design = get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + x2 +x3), sims)

test_that("Testing get design marix", {
  expect_equal(design$Y[, 1], sims$tstart, check.attributes = F, use.names = F)
  expect_equal(design$Y[, 2], sims$tstop, check.attributes = F, use.names = F)
  expect_equal(design$Y[, 3], sims$event, check.attributes = F, use.names = F)

  expect_equal(colnames(design$X), c("(Intercept)", "x1", "x2", "x3"))
  expect_equal(design$X[, 2:4], as.matrix(sims[, c("x1", "x2", "x3")]), check.attributes = F)
  })

########
# Test fixed terms

test_that("Fixed terms works as expected",{
  design_with_fixed <- get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3), sims)

  expect_equal(design_with_fixed$X[, c("x1", "x2", "x3")], as.matrix(sims[, c("x1", "x2", "x3")]),
               check.attributes = F)
  expect_length(design_with_fixed$fixed_terms, 0)


  design_with_fixed <- get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + ddFixed(x2) + x3), sims)
  expect_equal(design_with_fixed$X[, c("x1", "x3")], as.matrix(sims[, c("x1", "x3")]),
               check.attributes = F)
  expect_equal(as.matrix(design_with_fixed$fixed_terms), as.matrix(sims[, c("x2")]),
               check.attributes = F)


  design_with_fixed <- get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + ddFixed(x2) + ddFixed(x3)), sims)
  expect_equal(design_with_fixed$X[, c("x1"), drop = F], as.matrix(sims[, c("x1")]),
               check.attributes = F)
  expect_equal(as.matrix(design_with_fixed$fixed_terms), as.matrix(sims[, c("x2", "x3")]),
               check.attributes = F)

  dum <- function(x) cbind(x, -x)
  design_with_fixed <- get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + ddFixed(dum(x2)) + x3), sims)
  expect_equal(design_with_fixed$X[, c("x1", "x3")], as.matrix(sims[, c("x1", "x3")]),
               check.attributes = F)
  expect_equal(as.matrix(design_with_fixed$fixed_terms), dum(sims$x2),
               check.attributes = F)
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
    design = get_design_matrix(formula(survival::Surv(tstart, tstop, event) ~ x1 + x2 +x3), sims)
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
arg_list <- list(formula = formula(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3), data = sims, response = T)
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
arg_list <- list(formula = formula(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3), data = sims, response = T)
design <- do.call(get_design_matrix, arg_list)

arg_list$formula <- formula(survival::Surv(tstart, tstop, event) ~ -1 + x1 + x2 + x3)
design_no_intercept <- do.call(get_design_matrix, arg_list)

test_that("Design mat with vs. without intercept", {
  expect_equal(design$Y, design_no_intercept$Y)
  expect_equal(design$X[, 2:4], design_no_intercept$X[, 1:3])
  expect_equal(ncol(design_no_intercept$X[, 1:3]), 3)
})

test_that("is_for_discrete_model logic work",{
  by_ = 1
  t_max = 3

  data_ <- data.frame(matrix(
    dimnames = list(NULL, c("id", "tstart", "tstop", "event")),
    ncol = 4, byrow = T, data = c(
      # Id 1:
      #   is in bin 1 in either method
      #   has event
      1, 0, .5, 1,

      # Id 2:
      #   same as id 1
      2, 0, 1, 1,

      # Id 3:
      #   same as id 1 but with no event
      3, 0, .5, 0,

      # Id 4:
      #   same as id 2 but with no event
      4, 0, 1, 0,

      # Id 5:
      #   1 row in all three bins
      #   has event
      5, 0, 3, 1,

      # Id 6:
      #   5 row in all three bins
      #   no has event
      6, 0, 3, 0,

      # Id 7
      #   one row that is only in cont
      #   and is an event
      7, .25, .75, 1,

      # Id 8
      #   two rows where one cross two bins
      #   has event in the latter which is added two first in discrete
      8, 0, 1.25, 0,
      8, 1.25, 1.5, 1,

      # Id 9
      #   same as Id 7 but no event
      9, .25, .75, 0,

      # Id 10
      #   two rows where second cross two bins and has event after max_T
      10, .25, 1.5, 0,
      10, 1.5, 9, 1,

      # id 11
      #   same as 11 but with no event
      11, .25, 1.5, 0,
      11, 1.5, 9, 0,

      # id 12
      #   starts at time one and end at time 2 should be in second bin for both
      12, 1, 2, 0,

      # id 13
      #   same as id 12 but with an event
      13, 1, 2, 1,

      # id 14
      #   should be in both methods for last two bins. Has event in neither
      13, 1, 3.4, 1
  )))

  arg_list <- list(Y = survival::Surv(data_$tstart, data_$tstop, data_$event),
                   by = by_, max_T = t_max, id = data_$id)

  arg_list$is_for_discrete_model <- T
  get_risk_obj_discrete <- do.call(get_risk_obj, arg_list)

  expect_equal(get_risk_obj_discrete$is_event_in,
               c(0, 0, -1, -1, 2,
                 -1, -1, 1, -1, -1,
                 -1, -1, -1, -1, -1,
                 1, -1))

  expect_true(setequal(get_risk_obj_discrete$risk_sets[[1]], c(1:6, 8)))
  expect_true(setequal(get_risk_obj_discrete$risk_sets[[2]], c(5:6, 8, 11, 13, 15, 16, 17)))
  expect_true(setequal(get_risk_obj_discrete$risk_sets[[3]], c(5:6, 12, 14, 17)))

  arg_list$is_for_discrete_model <- F
  get_risk_obj_cont <- do.call(get_risk_obj, arg_list)

  expect_equal(get_risk_obj_cont$is_event_in,
               c(0, 0, -1, -1, 2,
                 -1, 0, -1, 1, -1,
                 -1, -1, -1, -1, -1,
                 1, -1))

  expect_true(setequal(get_risk_obj_cont$risk_sets[[1]], c(1:8, 10, 11, 13)))
  expect_true(setequal(get_risk_obj_cont$risk_sets[[2]], c(5:6, 8, 9, 11, 12, 13, 14, 15, 16, 17)))
  expect_true(setequal(get_risk_obj_cont$risk_sets[[3]], c(5:6, 12, 14, 17)))
})
