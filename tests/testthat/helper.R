suppressWarnings(RNGversion("3.5.0"))

# handle change of the default in all.equal
expect_equal <- function(...)
  testthat::expect_equal(..., check.environment = FALSE)

# Run before test
# See https://github.com/hadley/testthat/issues/544#issuecomment-260053774

# function to test that expression does not throw an error
expect_no_error = function(expr, env = parent.frame()){
  call <- match.call()
  eval(bquote(expect_error(.(call[[2]]), NA)), envir = env)
}

# expect_no_error(1 / "a")
# expect_no_error(1 / 1)

# load data example Fahrmeier (1994)
get_head_neck_cancer_data <- function(){
  is_censored = c(6, 27, 34, 36, 42, 46, 48:51,
                  51 + c(15, 30:28, 33, 35:37, 39, 40, 42:45))
  head_neck_cancer = data.frame(
    id = 1:96,
    start = rep(0, 96),
    stop = c(
      1, 2, 2, rep(3, 6), 4, 4, rep(5, 8),
      rep(6, 7), 7, 8, 8, 8,
      9, 9, 10, 10, 10, 11, 14, 14, 14, 15, 18, 18, 20, 20, 37,
      37, 38, 41, 45, 47, 47,

      2, 2, 3, rep(4, 4), rep(5, 5), rep(6, 5),
      7, 7, 7, 9, 10, 11, 12, 15, 16, 18, 18, 18, 21,
      21, 24, 25, 27, 36, 41, 44, 52, 54, 59, 59, 63, 67, 71, 76),
    event = !(1:96 %in% is_censored),
    group = factor(c(rep(1, 45 + 6), rep(2, 45))))

  head_neck_cancer$group = factor(head_neck_cancer$group, levels = c(2, 1))

  head_neck_cancer
}

########
# PBC data set from survival with the time-variying covariates
# See: https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf

get_pbc2_data <- function(){
  pbc <- survival::pbc
  pbcseq <- survival::pbcseq

  # avoid notes with CRAN tests
  id <- NULL
  sex <- NULL
  time <- NULL
  status <- NULL
  edema <- NULL
  age <- NULL
  event <- NULL
  tdc <- NULL
  day <- NULL
  albumin <- NULL
  protime <- NULL
  bili <- NULL

  temp <- subset(pbc, id <= 312, select=c(id, sex, time, status, edema, age))
  pbc2 <- survival::tmerge(
    temp, temp, id=id, death = event(time, status))
  pbc2 <- survival::tmerge(
    pbc2, pbcseq, id=id, albumin = tdc(day, albumin),
    protime = tdc(day, protime), bili = tdc(day, bili))
  pbc2 <- pbc2[, c("id", "tstart", "tstop", "death", "sex", "edema",
                   "age", "albumin", "protime", "bili")]

  pbc2
}

options(useFancyQuotes = FALSE)

if(interactive()){
  library(dynamichazard); library(testthat)

  list2env(as.list(asNamespace("dynamichazard")), envir = environment())

  if(!grepl("testtthat$", getwd()) && grepl("dynamichazard$", getwd()))
    setwd("tests/testthat")
}

# wrapper for expect_know_...
qu <- quote({
  file <- local({
    dir <- if(!interactive()) "./previous_results/" else
      paste0(gsub("(.+dynamichazard)(.*)", "\\1", getwd()),
             "/tests/testthat/previous_results/")

    paste0(dir, file)
  })})


do_update_tests <- FALSE
if(packageVersion("testthat") >= "1.0.2.9000"){
  trace(
    what = testthat::expect_known_value, at = 2, print = FALSE,
    tracer = qu)
  trace(
    what = testthat::expect_known_output, at = 2, print = FALSE,
    tracer = qu)

  expect_known_value <- testthat::expect_known_value
  formals(expect_known_value)$update <- do_update_tests
  expect_known_output <- testthat::expect_known_output
  formals(expect_known_output)$update <- do_update_tests

} else {
  trace(
    what = testthat::expect_output_file, at = 2, print = FALSE,
    tracer = qu)
  trace(
    what = testthat::expect_equal_to_reference, at = 2, print = FALSE,
    tracer = qu)

  expect_known_value <- testthat::expect_equal_to_reference
  formals(expect_known_value)$update <- do_update_tests
  expect_known_output <- testthat::expect_output_file
  formals(expect_known_output)$update <- do_update_tests

}

# And we may aswell just save it to a file
save_to_test <- function(obj, file_name, tolerance = sqrt(.Machine$double.eps)){
  if(!interactive())
    stop("save_to_test called not in interactive mode. Likely an error")

  cat("Largest sizes:\n")
  if(is.list(obj))
    print(head(sort(unlist(lapply(obj, object.size)), decreasing = T))) else
      print(object.size(obj))

  root <- gsub("(.+dynamichazard)(.*)", "\\1", getwd())
  out_file <- paste0(
    root, "/tests/testthat/previous_results/", file_name, ".RDS")
  saveRDS(obj, compress = T, out_file)

  cat("RDS file size is ", file.size(out_file) / 1000, "KB\n", sep = "")

  str_tol <- if(any_numeric(obj))
    paste0(", tolerance = ", signif(tolerance, 4)) else ""

  cat("Call 'expect_equal(", deparse(substitute(obj)), ", read_to_test(\"",  file_name, "\")",
      str_tol, ")' to test\n", sep = "")
}

any_numeric <- function(x){
  if(!is.recursive(x))
    return(is.numeric(x))

  for(i in x){
    out <- any_numeric(i)
    if(out)
      return(TRUE)
  }

  return(FALSE)
}

read_to_test <- function(file_name){
  readRDS(read_to_test_get_file_w_path(file_name))
}

read_to_test_get_file_w_path <- function(file_name){
  root <- gsub("(.+dynamichazard)(.*)", "\\1", getwd())
  path <- if(!interactive()) "./previous_results/" else
    paste0(root, "/tests/testthat/previous_results/")

  paste0(path, file_name, ".RDS")
}

test_if_file_exists <- function(file_name, test_expr){
  file_name <- gsub("(^.+)(\\.RDS)$", "\\1", file_name)

  if(!file.exists(file_w_path <- read_to_test_get_file_w_path(file_name))){
    cat("Skipped test as", sQuote(file_w_path), "does not exists\n")
  } else
    eval(substitute(test_expr), envir = parent.frame())
}

library(testthat)
library(dynamichazard)

options(ddhazard_max_threads = 2)
options(warn=1)

head_neck_cancer <- get_head_neck_cancer_data()
pbc2 <- get_pbc2_data()

if(exists("pbc", envir = environment(), inherits = FALSE))
  rm(pbc, envir = environment())
data("pbc", package = "survival")
pbc_org <- pbc
pbc <- pbc[, c("id", "time", "status", "age", "edema", "bili", "protime")]
pbc <- pbc[complete.cases(pbc), ]

# testbigglm_wrapper.R
library(biglm)
bigqr.init <- asNamespace("biglm")$bigqr.init

# Simulated data sets to test against
set.seed(6790753)

test_sim_func_logit <- dynamichazard:::test_sim_func_logit
get_sim <- function(n)
  test_sim_func_logit(n_series = n, n_vars = 10, t_0 = 0, t_max = 10,
                      x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(10),
                      is_fixed = 2:4,
                      intercept_start = -3, sds = c(.1, rep(.5, 10)))
logit_sim_200 <- get_sim(200)
logit_sim_500 <- get_sim(500)

# matplot(logit_sim_200$betas, type = "l", lty = 1)
# sum(logit_sim_200$res$event)
# hist(logit_sim_200$res$tstop[logit_sim_200$res$event == 1])

# matplot(logit_sim_500$betas, type = "l", lty = 1)
# sum(logit_sim_500$res$event)
# hist(logit_sim_500$res$tstop[logit_sim_500$res$event == 1])

set.seed(20406799)
test_sim_func_exp <- asNamespace("dynamichazard")$test_sim_func_exp

get_sim <- function(n)
  test_sim_func_exp(n_series = n, n_vars = 10, t_0 = 0, t_max = 10,
                      x_range = 1, x_mean = 0, re_draw = T, beta_start = rnorm(10),
                      is_fixed = 2:4,
                      intercept_start = -3, sds = c(.1, rep(.5, 10)))
exp_sim_200 <- get_sim(200)
exp_sim_500 <- get_sim(500)

# matplot(exp_sim_200$betas, type = "l", lty = 1)
# sum(exp_sim_200$res$event)
# hist(exp_sim_200$res$tstop[logit_sim_200$res$event == 1])

# matplot(exp_sim_500$betas, type = "l", lty = 1)
# sum(exp_sim_500$res$event)
# hist(exp_sim_500$res$tstop[logit_sim_500$res$event == 1])
