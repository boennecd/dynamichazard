# Run before test
# See https://github.com/hadley/testthat/issues/544#issuecomment-260053774

options(useFancyQuotes = FALSE)

if(interactive()){
  library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else source("R/test_utils.R")

  list2env(as.list(asNamespace("dynamichazard")), envir = environment())
}

# wrapper for expect_know_...
qu <- quote({
  file <- local({
    dir <- if(!interactive()) "./previous_results/" else
      paste0(stringr::str_match(getwd(), ".+dynamichazard"),
             "/tests/testthat/previous_results/")

    paste0(dir, file)
  })})

if(packageVersion("testthat") >= "1.0.2.9000"){
  trace(
    what = testthat::expect_known_value, at = 2, print = FALSE,
    tracer = qu)
  trace(
    what = testthat::expect_known_output, at = 2, print = FALSE,
    tracer = qu)

  expect_known_value <- testthat::expect_known_value
  formals(expect_known_value)$update <- FALSE
  expect_known_output <- testthat::expect_known_output
  formals(expect_known_output)$update <- FALSE

} else {
  trace(
    what = testthat::expect_output_file, at = 2, print = FALSE,
    tracer = qu)
  trace(
    what = testthat::expect_equal_to_reference, at = 2, print = FALSE,
    tracer = qu)

  expect_known_value <- testthat::expect_equal_to_reference
  formals(expect_known_value)$update <- FALSE
  expect_known_output <- testthat::expect_output_file
  formals(expect_known_output)$update <- FALSE

}

# And we may aswell just save it to a file
save_to_test <- function(obj, file_name, tolerance = sqrt(.Machine$double.eps)){
  if(!interactive())
    stop("save_to_test called not in interactive mode. Likely an error")

  cat("Largest sizes:\n")
  if(is.list(obj))
    print(head(sort(unlist(lapply(obj, object.size)), decreasing = T))) else
      print(object.size(obj))

  out_file <- paste0(
    stringr::str_match(getwd(), ".+dynamichazard"), "/tests/testthat/previous_results/", file_name, ".RDS")
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
  path <- if(!interactive()) "./previous_results/" else
    paste0(stringr::str_match(getwd(), ".+dynamichazard"), "/tests/testthat/previous_results/")

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

test_sim_func_logit <- asNamespace("dynamichazard")$test_sim_func_logit
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
