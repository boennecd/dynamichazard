context("Running test_cpp_utils.R:")

test_that("lambert_W0 gives close to same results as gsl::lambert_W0", {
  skip_if_not_installed("gsl")

  x <- log(seq(exp(1e-12), exp(3e-2), length.out = 1e3))

  expect_equal(sapply(x, lambert_W0_test), gsl::lambert_W0(x), tolerance = 5e-5)
})

test_that("trunc_lp_in_exponential_dist does not truncate when not needed", {
  skip_on_cran()
  .interactive <- interactive()

  eps <- exp(trunc_lp_in_exponential_dist_test_log_eps())

  test_func <- function(eta, .t, is_event){
    q <- bquote({
      ans <- trunc_lp_in_exponential_dist_test(
        eta = .(eta), at_risk_length = .(.t), is_event = .(is_event))
      expect_equal(
        ans$did_truncate,
        should_trunc <- .(eta) * .(is_event) - exp(.(eta)) * .(.t) < log(eps))

      if(.(.interactive)){
        cat(c(.(eta), .(.t), .(is_event)), "\n")
        cat("Likelihood is", exp(.(eta) * .(is_event) - exp(.(eta)) * .(.t)), "\n")
        cat("ans is ", unlist(ans), "\n")

      }
      if(!should_trunc){
        expect_equal(ans$exp_eta_trunc, exp(.(eta)))
        expect_equal(ans$eta_trunc, .(eta))

      } else{
        expect_equal(exp(ans$eta_trunc), ans$exp_eta_trunc)
        expect_equal(ans$eta_trunc * .(is_event) - exp(ans$eta_trunc) * .(.t), log(eps))

        if(.(is_event)){
          expect_equal(
            ans$eta_trunc < - exp(ans$eta_trunc) * .(.t),
            .(eta) < - exp(.(eta)) * .(.t))
        }

      }
    })

    eval(q, envir = environment())
  }

  test_vals <- expand.grid(
    eta = -25:25,
    t = c(1e-2, 1, 1e2),
    is_event = c(TRUE, FALSE))

  invisible(
    mapply(test_func, eta = test_vals$eta, .t = test_vals$t, is_event = test_vals$is_event)
  )
})


test_that("round_if_almost_eq rounds to nearest boundary as expected", {
  round_if_almost_eq <- asNamespace("dynamichazard")$round_if_almost_eq
  set.seed(56219385)
  n <- 1e2

  x <- replicate(n, {
    x <- 1
    for(i in 1:10){
      tmp <- rexp(1, 1)
      n <- sample.int(1e2, 1)
      delta <- tmp / n
      for(j in 1:n)
      x <- x - delta
      x <- x + tmp
    }
    x
  })

  boundaries <- seq(0, 2, .1)
  x_ord <- order(x) - 1L
  x_trans <- drop(round_if_almost_eq(x, x_ord, boundaries = boundaries))

  correct_bound <- boundaries[11]
  expect_gt(sx <- sum(x_trans == correct_bound), sum(x == correct_bound))
  expect_equal(sx, length(x))

  # values not at the boundary should not get truncated
  x <- rnorm(n)
  x_ord <- order(x) - 1L
  x_trans <- drop(round_if_almost_eq(x, x_ord, boundaries = boundaries))
  expect_equal(x_trans, x)
})
