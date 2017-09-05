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


