context("Testing design_mat_and_risk_obj")

# Simulate data
sims = logit_sim_200$res

#######
# Test design mat
design = get_design_matrix(
  survival::Surv(tstart, tstop, event) ~ x1 + x2 +x3, sims)

test_that("Testing get design marix", {
  expect_equal(design$Y[, 1], sims$tstart, check.attributes = F, use.names = F)
  expect_equal(design$Y[, 2], sims$tstop, check.attributes = F, use.names = F)
  expect_equal(design$Y[, 3], sims$event, check.attributes = F, use.names = F)

  expect_equal(colnames(design$X), c("(Intercept)", "x1", "x2", "x3"))
  expect_equal(design$X[, 2:4], as.matrix(sims[, c("x1", "x2", "x3")]), check.attributes = F)
  })

########
# Test fixed terms

test_that("different ways of fixing the intercept gives the same result", {
  .data <- data.frame(
    y = rnorm(25) > 0, x = sample(letters[1:3], 25, replace = TRUE),
    tstop = rexp(25, 1))
  test_formulas <- list(
    Surv(tstop, y) ~ ddFixed_intercept() + x,
    Surv(tstop, y) ~ -1 + ddFixed_intercept() + x)

  test_exp <- expression({
    results <- lapply(test_formulas, get_design_matrix, data = .data)

    for(i in 2:length(test_formulas))
      eval(substitute({
        expect_equal(results[[1]], results[[i]])
        r_new <- with(results[[i]], get_design_matrix(
          formula = test_formulas[[i]],
          data = .data[1:10, ], Terms = terms, xlev = xlev,
          has_fixed_intercept = has_fixed_intercept))
        expect_equal(results[[i]][["X"]][1:10, , drop = FALSE], r_new[["X"]])
        expect_equal(results[[i]][["fixed_terms"]][1:10, , drop = FALSE],
                     r_new[["fixed_terms"]])
        expect_equal(results[[i]][["Y"]][1:10, ], r_new[["Y"]][1:10, ])
      }, list(i = i)))
  })
  eval(test_exp)

  # Key thing about poly is that forms an orthogonal basis so we need to use
  # predvars in the terms object to get the same resutls (if I am correct)
  .data$x <- rnorm(25)
  test_formulas <- list(
    Surv(tstop, y) ~ ddFixed_intercept() + poly(x, degree = 3),
    Surv(tstop, y) ~ -1 + ddFixed_intercept() + poly(x, degree = 3))
  eval(test_exp)

  test_formulas <- list(
    Surv(tstop, y) ~ ddFixed_intercept() + x,
    Surv(tstop, y) ~ -1 + ddFixed_intercept() + x)
  eval(test_exp)
})

test_that(
  "Previous way of fixing the intercept gives the correct error message", {
    .data <- data.frame(tstop = rep(1, 10), event = c(rep(0, 5), rep(1, 5)),
                        x = 1:10)

    lapply(list(
      Surv(tstop, event) ~ ddFixed(1) + x1,
      Surv(tstop, event) ~ -1 + ddFixed(1) + x1,
      Surv(tstop, event) ~ ddFixed(rep(1, nrow(.data)))),
      function(frm) eval(bquote(
        expect_error(
          get_design_matrix((frm), .data),
          "^All elements in call to .ddFixed. is one. This is potentially a"
          ))))
  })

test_that("Fixed terms works as expected",{
  dum <- function(x) cbind(x, -x)

  design_with_fixed <- get_design_matrix(
    formula(survival::Surv(tstart, tstop, event) ~
              x1 + ddFixed(dum(x2)) + x3), sims)

  expect_equal(design_with_fixed$X[, c("x1", "x3")], as.matrix(sims[, c("x1", "x3")]),
               check.attributes = F)
  expect_equal(as.matrix(design_with_fixed$fixed_terms), dum(sims$x2),
               check.attributes = F)
})

test_that("Fixed terms works as expected",{
  design_with_fixed <- get_design_matrix(
    formula(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3), sims)

  expect_equal(
    design_with_fixed$X[, c("x1", "x2", "x3")],
    as.matrix(sims[, c("x1", "x2", "x3")]), check.attributes = F)
  expect_length(design_with_fixed$fixed_terms, 0)


  design_with_fixed <- get_design_matrix(formula(
    survival::Surv(tstart, tstop, event) ~ x1 + ddFixed(x2) + x3), sims)
  expect_equal(design_with_fixed$X[, c("x1", "x3")],
               as.matrix(sims[, c("x1", "x3")]),
               check.attributes = F)
  expect_equal(as.matrix(design_with_fixed$fixed_terms),
               as.matrix(sims[, c("x2")]),
               check.attributes = F)


  design_with_fixed <- get_design_matrix(
    survival::Surv(tstart, tstop, event) ~ x1 + ddFixed(x2) + ddFixed(x3),
    sims)
  expect_equal(design_with_fixed$X[, c("x1"), drop = F],
               as.matrix(sims[, c("x1")]),
               check.attributes = F)
  expect_equal(as.matrix(design_with_fixed$fixed_terms),
               as.matrix(sims[, c("x2", "x3")]),
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

        expect_true(all(
          (design$Y[not_set, 1] > t_ - by_ |     # after interval start
             design$Y[not_set, 2] <= t_ - by_) | # stop before interval start
            (design$Y[not_set, 1] <= t_ - by_ &  # is inside a bin
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
    if(sum(inside_of_bin_and_event) > 0)
      test_that("All those inside of bins with events have has the correct previous row as event", {
        ids <- sims$id[inside_of_bin_and_event]
        ids_to_check <- sample(ids, size = min(1e3, length(ids)), replace = F)

        for(id in ids_to_check){
          rows_ <- sims[ids_to_row_indicies[[id]], ]
          event_time <- max(rows_$tstop)
          bin_start <- min(
            risk_set$event_times[risk_set$event_times >= event_time]) - by_

          should_be_event <- which(
            rows_$tstart <= bin_start &  bin_start < rows_$tstop)
          expect_true(
            risk_set$is_event_in[
              ids_to_row_indicies[[id]][should_be_event]] >= 0)
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
    if(sum(inside_of_bin_and_event) > 0)
      test_that("All those inside of bins with events have has the correct previous row as event", {
        ids <- sims$id[inside_of_bin_and_event]
        ids_to_check <- sample(ids, size = min(1e3, length(ids)), replace = F)
        for(id in ids_to_check){
          rows_ <- sims[ids_to_row_indicies[[id]], ]
          event_time <- max(rows_$tstop)
          bin_start <- min(
            risk_set$event_times[risk_set$event_times >= event_time]) - by_

          should_be_event <- which(
            rows_$tstart <= bin_start &  bin_start < rows_$tstop)

          expect_true(risk_set$is_event_in[
            ids_to_row_indicies[[id]][should_be_event]] >= 0)
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
  expect_null(design_no_Y$Y)
})

arg_list$data <- sims[, !colnames(sims) %in% c("tstart", "tstop", "event")]

tryCatch({
  test_that("calling get_design_matrix without response and also no response in data works and give the same as calling with response",{
    design_no_Y_also_in_data <- do.call(get_design_matrix, arg_list)

    expect_equal(design_no_Y_also_in_data$X, design_no_Y$X)
    expect_equal(design_no_Y_also_in_data$formula, design_no_Y$formula)
    expect_null(design_no_Y_also_in_data$Y)
  })
}, error = function(...){
  test_that("calling get_design_matrix without response and also no response in data",
            expect_true(FALSE))
})

test_that("dimension are correct with get_design_matrix with different combinations of fixed or not fixed factor",{
  #####
  # Only time varying w/ no intercept
  dsng_mat <- get_design_matrix(
    Surv(tstart, tstop, event) ~
      -1 + as.factor(x1 > 0) + as.factor(x2 > 0) + as.factor(x3 > 0), sims)

  expect_equal(ncol(dsng_mat$fixed_terms), 0)
  expect_equal(ncol(dsng_mat$X), 4)
  expect_equal(colnames(dsng_mat$X),
               c("as.factor(x1 > 0)FALSE", "as.factor(x1 > 0)TRUE",
                 "as.factor(x2 > 0)TRUE", "as.factor(x3 > 0)TRUE"))

  #####
  # Only time varying w/ intercept
  dsng_mat <- get_design_matrix(
    Surv(tstart, tstop, event) ~
      as.factor(x1 > 0) + as.factor(x2 > 0) + as.factor(x3 > 0), sims)

  expect_equal(ncol(dsng_mat$fixed_terms), 0)
  expect_equal(ncol(dsng_mat$X), 4)
  expect_equal(colnames(dsng_mat$X),
               c("(Intercept)", "as.factor(x1 > 0)TRUE",
                 "as.factor(x2 > 0)TRUE", "as.factor(x3 > 0)TRUE"))

  #####
  # No intercept w/ mixed fixed and time-varying effects
  dsng_mat <- get_design_matrix(
    Surv(tstart, tstop, event) ~
      -1 + ddFixed(as.factor(x1 > 0)) + as.factor(x2 > 0) +
      as.factor(x3 > 0), sims)

  expect_equal(ncol(dsng_mat$fixed_terms), 2)
  expect_equal(colnames(dsng_mat$fixed_terms),
               c("ddFixed(as.factor(x1 > 0))FALSE", "ddFixed(as.factor(x1 > 0))TRUE"))
  expect_equal(ncol(dsng_mat$X), 2)
  expect_equal(colnames(dsng_mat$X),
               c("as.factor(x2 > 0)TRUE", "as.factor(x3 > 0)TRUE"))

  #####
  # Fixed intercept w/ mixed fixed and time-varying effects
  dsng_mat <- get_design_matrix(
    Surv(tstart, tstop, event) ~
      ddFixed_intercept() + ddFixed(as.factor(x1 > 0)) + as.factor(x2 > 0) +
      as.factor(x3 > 0), sims)
  expect_equal(ncol(dsng_mat$fixed_terms), 2)
  expect_equal(colnames(dsng_mat$fixed_terms),
               c("(Intercept)", "ddFixed(as.factor(x1 > 0))TRUE"))
  expect_equal(ncol(dsng_mat$X), 2)
  expect_equal(colnames(dsng_mat$X),
               c("as.factor(x2 > 0)TRUE", "as.factor(x3 > 0)TRUE"))
})

# Test that removal of intercepts works as in lm
arg_list <- list(formula = survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3,
                 data = sims, response = T)
design <- do.call(get_design_matrix, arg_list)

arg_list$formula <- survival::Surv(tstart, tstop, event) ~ -1 + x1 + x2 + x3
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
      #   same as id 1 but with no event. Hence not in
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

  expect_true(setequal(get_risk_obj_discrete$risk_sets[[1]], c(c(1, 2, 4, 5, 6), 8)))
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

#####
# "By hand" test of risk_obj

# This is a useful illustration (I find)
plot_exp <- expression({
  id_unique <- unique(id)

  ylim <- c(1, length(id_unique))
  xlim <- range(tstart, tstop)
  plot(xlim, ylim + c(-.5, .5), type = "n", ylab = "", xlab = "time")
  abline(v = 0:ceiling(xlim[2]), lty = 2)

  for(y in 1:ylim[2]){
    i <- id_unique[y]
    is_in <- which(id == i)
    n <- length(is_in)

    tsta <- tstart[is_in]
    tsto <- tstop[is_in]

    pch <- ifelse(event[is_in] == 1, rep(16,  n), rep(4, n))

    for(u in 1:n){
      lines(c(tsta[u], tsto[u]), rep(y, 2))
    }
    points(tsto, rep(y, n), pch = pch)
    text(tsta - .05, rep(y + .5, n), labels = is_in)
  }
})

# The test data
test_dat <- data.frame(matrix(
  byrow = T, ncol = 4,
  dimnames = list(NULL, c("tstart", "tstop", "event", "id")), c(
    0,   1, 0, 1,
    1,   2, 1, 1, # 2
    0,   1, 0, 2,
    1,   2, 0, 2,

    0,   2, 0, 3,
    0,   2, 1, 4, # 6

    1,   2, 0, 5,
    1,   2, 1, 6, # 8

    0, 1.5, 1, 7, # 9
    0, 1.5, 0, 8,

    0,  .5, 0,  9,
    0,  .5, 1, 10, # 12

    .25, .5, 0, 11,
    .25, .5, 1, 12,

     .5, 1.3, 0, 13,
    1.3, 1.5, 0, 13,
     .5, 1.3, 0, 14, # 17
    1.3, 1.5, 1, 14,

      0, 1.1, 0, 15,
    1.3, 1.6, 1, 15,
      0, 1.1, 0, 16,
    1.3, 1.6, 0, 16
  )))

## Comment back to get plot
# with(test_dat, eval(plot_exp))

risk_obj <- with(test_dat, get_risk_obj(
  Y = Surv(tstart, tstop, event),
  by = 1, max_T = 2, id = id, is_for_discrete_model = T))

test_that("Each risk set contain the expect indvidual with the correct rows", {
  risk_sets_by_ids <- lapply(risk_obj$risk_sets, function(x) test_dat$id[x])
  expect_true(setequal(risk_sets_by_ids[[1]], c(1:4, 7, 8, 10, 15, 16)))
  expect_true(setequal(sort(risk_sets_by_ids[[2]]), c(1:7, 14, 15)))

  expect_equal(sort(risk_obj$risk_sets[[1]]), c(1, 3, 5, 6, 9, 10, 12, 19, 21))
  expect_equal(sort(risk_obj$risk_sets[[2]]), c(2, 4, 5, 6, 7, 8, 9, 17, 19))
})

test_that("Individuals are markeds as events in the correct bins", {
  is_ev <- c(2, 6, 8, 9, 12, 17, 19)

  expect_true(all(risk_obj$is_event_in[-is_ev] == -1))
  expect_equal(risk_obj$is_event_in[is_ev], c(1, 1, 1, 1, 0, 1, 1))
})

#####
# Test randomize order

set.seed(875080)
X_Y <- get_design_matrix(
  Surv(tstart, tstop, death == 2) ~ ddFixed(age) + albumin, pbc2)
risk_obj_org <- risk_obj <- get_risk_obj(X_Y$Y, 100, max_T = 3600, pbc2$id)
X_Y[1:2] <- lapply(X_Y[1:2], data.frame)
X_Y_org <- X_Y

w_org <- w <- sample(1:3, replace = T, nrow(pbc2))

eval(get_permu_data_exp(X_Y[1:3], risk_obj, w))


test_that("Permutating gives the data frame and risk set permutated", {
  got_plyr <- requireNamespace("plyr", quietly = TRUE)

  for(i in 1:2){
    expect_true(any(as.matrix(X_Y[[i]]) != as.matrix(X_Y_org[[i]])))
    expect_equal(nrow(X_Y[[i]]), nrow(X_Y_org[[i]]))
    if(got_plyr)
      expect_equal(nrow(suppressMessages(plyr::match_df(X_Y_org[[i]], X_Y[[i]]))),
        nrow(X_Y_org[[i]]))
  }
  expect_true(sum(is.na(match(X_Y_org$Y, X_Y$Y))) == 0)

  expect_true(setequal(w_org, w))
  expect_true(any(w_org != w))

  expect_equal(length(risk_obj$risk_sets), length(risk_obj_org$risk_sets))
  for(i in seq_along(risk_obj$risk_sets)){
    r1 <- risk_obj$risk_sets[[i]]
    r2 <- risk_obj_org$risk_sets[[i]]

    expect_true(!is.unsorted(risk_obj$risk_sets[[i]]))
    expect_true(any(r1 != r2))
    expect_equal(length(r1), length(r2))
  }

  expect_true(any(risk_obj$is_event_in != risk_obj_org$is_event_in))
  expect_true(setequal(risk_obj$is_event_in, risk_obj_org$is_event_in))
})

eval(get_permu_data_rev_exp(X_Y[1:3], risk_obj, w))

test_that("Reverse permutating the data frame and the initial values", {
  for(i in 1:3)
    expect_equal(as.matrix(X_Y[[i]]), as.matrix(X_Y_org[[i]]))

  for(i in seq_along(risk_obj$risk_sets)){
    r1 <- risk_obj$risk_sets[[i]]
    r2 <- risk_obj_org$risk_sets[[i]]

    expect_equal(r1, r2)
  }

  expect_equal(w, w_org)

  expect_equal(risk_obj$is_event_in, risk_obj_org$is_event_in)
})

#####
# Test sorting
set.seed(875080 + 1)
X_Y <- get_design_matrix(
  Surv(tstart, tstop, death == 2) ~ ddFixed(age) + albumin, pbc2)
risk_obj_org <- risk_obj <- get_risk_obj(X_Y$Y, 100, max_T = 3600, pbc2$id)
X_Y[1:2] <- lapply(X_Y[1:2], data.frame)
X_Y_org <- X_Y

w_org <- w <- sample(1:3, replace = T, nrow(pbc2))

eval(get_order_data_exp(X_Y[1:3], risk_obj, w))

test_that("Ordering gives the data frame and risk set orded", {
  got_plyr <- requireNamespace("plyr", quietly = TRUE)

  for(i in 1:2){
    expect_true(any(as.matrix(X_Y[[i]]) != as.matrix(X_Y_org[[i]])))
    expect_equal(nrow(X_Y[[i]]), nrow(X_Y_org[[i]]))
    if(got_plyr){
      expect_equal(nrow(suppressMessages(plyr::match_df(X_Y_org[[i]], X_Y[[i]]))),
                   nrow(X_Y_org[[i]]))
    }
  }
  expect_true(sum(is.na(match(X_Y_org$Y, X_Y$Y))) == 0)

  expect_true(setequal(w_org, w))
  expect_true(any(w_org != w))

  expect_equal(length(risk_obj$risk_sets), length(risk_obj_org$risk_sets))
  for(i in seq_along(risk_obj$risk_sets)){
    r1 <- risk_obj$risk_sets[[i]]
    r2 <- risk_obj_org$risk_sets[[i]]

    expect_true(!is.unsorted(risk_obj$risk_sets[[i]]))
    expect_true(any(r1 != r2))
    expect_equal(length(r1), length(r2))
  }

  expect_true(any(risk_obj$is_event_in != risk_obj_org$is_event_in))
  expect_true(setequal(risk_obj$is_event_in, risk_obj_org$is_event_in))

  expect_false(is.unsorted(X_Y$Y[, 2]))
})

eval(get_order_data_rev_exp(X_Y[1:3], risk_obj, w))

test_that("Reverse permutating the data frame and the initial values", {
  for(i in 1:3)
    expect_equal(as.matrix(X_Y[[i]]), as.matrix(X_Y_org[[i]]))

  for(i in seq_along(risk_obj$risk_sets)){
    r1 <- risk_obj$risk_sets[[i]]
    r2 <- risk_obj_org$risk_sets[[i]]

    expect_equal(r1, r2)
  }

  expect_equal(w, w_org)

  expect_equal(risk_obj$is_event_in, risk_obj_org$is_event_in)

  expect_true(is.unsorted(X_Y$Y[, 2]))
})

#####
# By hand test of risk_obj with parallal

test_that("Parallel and non-parallel version gives the same for get_risk_set. The former gives a warning", {
  set.seed(4296745)
  s <-
    test_sim_func_logit(
      n_series = 5000, n_vars = 1, beta_start = 1,
      intercept_start = - 5, sds = c(sqrt(.1), rep(1, 1)),
      x_range = 1, x_mean = .5)$res

  p1 <- with(s, get_risk_obj(
    Y = Surv(tstart, tstop, event),
    by = 1, max_T = 10, id = id, is_for_discrete_model = T))

  expect_warning(
    p2 <- with(s, get_risk_obj(
      Y = Surv(tstart, tstop, event),
      by = 1, max_T = 10, id = id, is_for_discrete_model = T,
      n_threads = 2, min_chunk = 1000)),
    "'n_threads' greater than one is no longer supported", fixed = TRUE)

  for(i in 1:10)
    expect_true(setequal(p1$risk_sets[[i]],
                         p2$risk_sets[[i]]))

  expect_equal(p1[-1], p2[-1])
})
