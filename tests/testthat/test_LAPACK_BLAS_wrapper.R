if(interactive()){
  chol_rank_one_update <- with(environment(ddhazard), chol_rank_one_update)

  square_tri_inv <- with(environment(ddhazard), square_tri_inv)

  symmetric_mat_chol <- with(environment(ddhazard), symmetric_mat_chol)

  tri_mat_times_vec <- with(environment(ddhazard), tri_mat_times_vec)
}

test_name <- "Testing LAPACK wrapper functions"
cat("\nRunning", test_name, "\n")

set.seed(3144)

# Function to simulate the a covariance matrix
get_random_sym_post_def_mat <- function(n){
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp)
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  t(O) %*% diag(runif(n, 0, 10)) %*% O
}


#####
# Test for rank one update

test_that("rank one update of chol decomp works", {
  for(n in c(5, 10, 50, 100)){
    r_mat <- get_random_sym_post_def_mat(n)

    d1 <- t(chol(r_mat))
    x <- rnorm(n)
    d2 <- t(chol(r_mat + x %*% t(x)))

    d1_tmp <- d1 # Take copy as it will be overwritten

    chol_rank_one_update(d1_tmp, x)

    expect_equal(d2, d1_tmp)
  }
})

# # Should be rougly quadratic in n
# (ns <- 2^(5:11))
# in_dat <- lapply(
#   ns, function(x) get_random_sym_post_def_mat(x))
# xs <- lapply(ns, rnorm)
#
# library(microbenchmark)
# out <- microbenchmark(
#   n32 = {
#     x <- xs[[1]]
#     tmp <- in_dat[[1]]
#     chol_rank_one_update(tmp, x)
#   },
#
#   n64 = {
#     x <- xs[[2]]
#     tmp <- in_dat[[2]]
#     chol_rank_one_update(tmp, x)
#   },
#
#   n128 = {
#     x <- xs[[3]]
#     tmp <- in_dat[[3]]
#     chol_rank_one_update(tmp, x)
#   },
#
#   n256 = {
#     x <- xs[[4]]
#     tmp <- in_dat[[4]]
#     chol_rank_one_update(tmp, x)
#   },
#
#   n512 = {
#     x <- xs[[5]]
#     tmp <- in_dat[[5]]
#     chol_rank_one_update(tmp, x)
#   },
#
#   n1024 = {
#     x <- xs[[6]]
#     tmp <- in_dat[[6]]
#     chol_rank_one_update(tmp, x)
#   },
#
#   n2048 = {
#     x <- xs[[7]]
#     tmp <- in_dat[[7]]
#     chol_rank_one_update(tmp, x)
#   })
#
# out_sum <- summary(out)
# plot(out_sum$median ~ ns , log="xy")
# lm(log(out_sum$median) ~ log(ns))$coefficients
# lm(log(out_sum$median) ~ log(ns), subset = 3:length(ns))$coefficients

#####
# Test for square tri inv

test_that("Square triangular inversion works followed by rank one update", {
  for(n in c(5, 10, 50, 100)){
    r_mat <- get_random_sym_post_def_mat(n)

    d1 <- t(chol(r_mat))
    d2 <- matrix(0, ncol = n, nrow = n)

    square_tri_inv(d1, d2)

    expect_equal(solve(d1), d2)
  }
})

#####
# Test of symmetric matrix Choleksy decomposition

test_that("Symmetric matrix Choleksy decomposition works", {
  for(n in c(5, 10, 50, 100)){
    r_mat <- get_random_sym_post_def_mat(n)

    d1 <- t(chol(r_mat))
    d2 <- matrix(0, ncol = n, nrow = n)

    symmetric_mat_chol(r_mat, d2)

    expect_equal(d1, d2)
  }
})

#####

test_that("Triangular matrix times vector works", {
  for(n in c(10, 50, 100)){
    r_mat <- get_random_sym_post_def_mat(n)
    d1 <- t(chol(r_mat))

    v1 <- rnorm(n)
    out <- rep(0, n)
    tri_mat_times_vec(d1, v1, out, F)
    expect_equal(c(d1 %*% v1), out)

    v2 <- rnorm(n)
    out <- rep(0, n)
    tri_mat_times_vec(d1, v2, out, T)
    expect_equal(c(t(d1) %*% v2), out)

    v3 <- rnorm(n / 2)
    out <- rep(0, n)
    tri_mat_times_vec(d1, v3, out, F)
    expect_equal(c(d1[, 1:(n/2)] %*% v3), out)

    v4 <- rnorm(n / 2)
    out <- rep(0, n)
    tri_mat_times_vec(d1, v4, out, T)
    expect_equal(c(t(d1)[, 1:(n/2)] %*% v4), out)
  }
})


# Had issues with win builder. Thus, these lines
cat("\nFinished", test_name, "\n")
