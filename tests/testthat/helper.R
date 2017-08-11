# Run before test
# See https://github.com/hadley/testthat/issues/544#issuecomment-260053774

if(interactive()){
  library(dynamichazard); library(testthat)

  if(grepl("testthat$", getwd()))
    source("../../R/test_utils.R") else
      source("R/test_utils.R")

  exp_model_names <- with(environment(ddhazard), exp_model_names)

  # test_LAPACK_BLAS_wrapper.R
  chol_rank_one_update <- with(environment(ddhazard), chol_rank_one_update)
  square_tri_inv <- with(environment(ddhazard), square_tri_inv)
  symmetric_mat_chol <- with(environment(ddhazard), symmetric_mat_chol)
  tri_mat_times_vec <- with(environment(ddhazard), tri_mat_times_vec)
  sym_mat_rank_one_update <- with(environment(ddhazard), sym_mat_rank_one_update)
  solve_w_precomputed_chol_test <- with(environment(ddhazard), solve_w_precomputed_chol_test)

  # test_SMA.R
  SMA_hepler_logit_compute_length <-
    with(environment(ddhazard), SMA_hepler_logit_compute_length)
  SMA_hepler_logit_second_d <-
    with(environment(ddhazard), SMA_hepler_logit_second_d)
  SMA_hepler_exp_compute_length <-
    with(environment(ddhazard), SMA_hepler_exp_compute_length)
  SMA_hepler_exp_second_d <-
    with(environment(ddhazard), SMA_hepler_exp_second_d)

  # testbigglm_wrapper.R
  bigglm_updateQR_rcpp <- with(environment(ddhazard), bigglm_updateQR_rcpp)
  bigglm_regcf_rcpp <- with(environment(ddhazard), bigglm_regcf_rcpp)

  # testboot_est.R
  get_frac_n_weights <- with(environment(ddhazard), get_frac_n_weights)

  # testdesign_mat_and_risk_obj.R
  get_permu_data_exp <-  with(environment(ddhazard), get_permu_data_exp)
  get_permu_data_rev_exp <- with(environment(ddhazard), get_permu_data_rev_exp)
  get_order_data_exp <-  with(environment(ddhazard), get_order_data_exp)
  get_order_data_rev_exp <- with(environment(ddhazard), get_order_data_rev_exp)
  get_design_matrix <- environment(ddhazard)$get_design_matrix

  # teststatic_glm.R
  get_design_matrix <- environment(ddhazard)$get_design_matrix
}

# hint use this regexp '\r\n\s*\[\d+\]\s' and replace with '\r\n,' followed by
# as search for '(?<=\d)\s+(?=\d)' and replace with ','
# or be lazy and use this function
str_func <- function(x, n_digist = 16){
  tmp <- capture.output(print(unname(c(x)), digits = n_digist))
  tmp <- sapply(tmp, gsub, pattern = "\\s*\\[1\\]\\s", replacement = "", USE.NAMES = F)
  tmp <- sapply(tmp, gsub, pattern = "\\s*\\[\\d+\\]\\s", replacement = ",\\ ", USE.NAMES = F)
  tmp <- sapply(tmp, gsub, pattern = "(?<=(\\d)|(NA))\\s+(?=(\\d)|-|(NA))", replacement = ",\\ ", perl = T, USE.NAMES = F)

  tmp <- paste0(c("c(", tmp, " )"), collapse = "")

  max_lengt <- floor(8191 * .75)
  n_nums_before_break = floor(max_lengt / (n_digist + 4))
  gsub(pattern = paste0("((\\d,\\ .*?){", n_nums_before_break - 1, "})(,\\ )"), replacement = "\\1,\n\\ ",
       x = tmp, perl = T)
}

get_expect_equal <- function(x, eps, file = ""){
  op_old <- options()
  on.exit(options(op_old))
  options(max.print = 1e7)

  arg_name <- deparse(substitute(x))
  expects <- unlist(lapply(x, str_func))
  tol_string = if(!missing(eps)) paste0("\n, tolerance = " , eval(bquote(.(eps)))) else ""
  expects <- mapply(function(e, index_name)
    paste0("expect_equal(unname(c(", arg_name, "$", index_name, ")),\n", e,
           tol_string, ")", collapse = ""),
    e = expects, index_name = names(expects))

  out <- paste0(c("{", paste0(expects, collapse = "\n\n"), "}\n"), collapse = "\n")
  cat(out, file = file)
  invisible()
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
  path <- if(!interactive()) "./previous_results/" else
    paste0(stringr::str_match(getwd(), ".+dynamichazard"), "/tests/testthat/previous_results/")

  readRDS(paste0(path, file_name, ".RDS"))
}

library(testthat)
library(dynamichazard)

options(ddhazard_max_threads = 2)
options(warn=1)

head_neck_cancer <- get_head_neck_cancer_data()
pbc2 <- get_pbc2_data()

data("pbc", package = "survival")
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

# ######
# # Debugging what takes time
# files <- list.files("tests/testthat")
# files <- files[grepl("^test", files)]
#
# time_taken <- sapply(files, function(f){
#   cat("Running", sQuote(f), "\n")
#   print(out <- system.time(testthat::test_file(paste0("tests/testthat/", f))))
#   cat("\n")
#
#   out
# })
#
# time_taken <- t(time_taken)
# time_taken[order(time_taken[, "elapsed"], decreasing = TRUE), ]
# sum(time_taken[, "elapsed"])

# ######
# # To test a particular file
# system.time(testthat::test_file("tests/testthat/testpredict.R"))
#
# Or use:
# test_that <- function(desc, code){
#   cat("\nRunning", sQuote(desc), "\n")
#   .time <- system.time(out <- testthat::test_that(desc, code))
#   print(.time)
#   cat("\n")
#
#   out
# }
