context("Testing IWLS")

# Use the example here
test_data <- if(grepl("testthat", getwd()))
  read.csv("binary.csv") else
    read.csv("tests/testthat/binary.csv")

# GLM fit
form <- formula(admit ~ .)
X <- model.matrix(form, data = test_data)
y <- model.response(model.frame(form, data = test_data))
beta_0 <- rep(0, ncol(X))

# Compare fits
test_that("glm.fit and IWLS return the same for binomial model", {
  for(is_random_offset in c(F, T)){
    offsets <- if(is_random_offset)
      rexp(nrow(test_data), 1) else rep(0, nrow(test_data))

    glm_fit <- glm.fit(x = X, y = y, start = beta_0, offset = offsets, family = binomial())

    IWLS_fit <-
      IWLS_logit(X = X, y = y, beta = rep(0, ncol(X)), offsets = offsets, it_max = 1e2,
                 eps = glm.control()$epsilon)


    expect_equal(unname(glm_fit$coefficients), c(IWLS_fit))
}})

# Do the same for poisson models
# data from http://www.theanalysisfactor.com/generalized-linear-models-in-r-part-6-poisson-regression-count-variables/
test_data <-
  structure(list(Days = c(1L, 2L, 3L, 3L, 4L, 4L, 4L, 6L, 7L, 8L,
                          8L, 8L, 8L, 12L, 14L, 15L, 17L, 17L, 17L, 18L, 19L, 19L, 20L,
                          23L, 23L, 23L, 24L, 24L, 25L, 26L, 27L, 28L, 29L, 34L, 36L, 36L,
                          42L, 42L, 43L, 43L, 44L, 44L, 44L, 44L, 45L, 46L, 48L, 48L, 49L,
                          49L, 53L, 53L, 53L, 54L, 55L, 56L, 56L, 58L, 60L, 63L, 65L, 67L,
                          67L, 68L, 71L, 71L, 72L, 72L, 72L, 73L, 74L, 74L, 74L, 75L, 75L,
                          80L, 81L, 81L, 81L, 81L, 88L, 88L, 90L, 93L, 93L, 94L, 95L, 95L,
                          95L, 96L, 96L, 97L, 98L, 100L, 101L, 102L, 103L, 104L, 105L,
                          106L, 107L, 108L, 109L, 110L, 111L, 112L, 113L, 114L, 115L),
                 Students = c(6L, 8L, 12L, 9L, 3L, 3L, 11L, 5L, 7L, 3L, 8L,
                              4L, 6L, 8L, 3L, 6L, 3L, 2L, 2L, 6L, 3L, 7L, 7L, 2L, 2L, 8L,
                              3L, 6L, 5L, 7L, 6L, 4L, 4L, 3L, 3L, 5L, 3L, 3L, 3L, 5L, 3L,
                              5L, 6L, 3L, 3L, 3L, 3L, 2L, 3L, 1L, 3L, 3L, 5L, 4L, 4L, 3L,
                              5L, 4L, 3L, 5L, 3L, 4L, 2L, 3L, 3L, 1L, 3L, 2L, 5L, 4L, 3L,
                              0L, 3L, 3L, 4L, 0L, 3L, 3L, 4L, 0L, 2L, 2L, 1L, 1L, 2L, 0L,
                              2L, 1L, 1L, 0L, 0L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 0L, 0L,
                              0L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)), .Names = c("Days", "Students"
                              ), class = "data.frame", row.names = c(NA, -109L))

IWLS_poisson <- with(environment(ddhazard), IWLS_poisson)

form <- formula(Students ~ Days)
X <- model.matrix(form, data = test_data)
y <- model.response(model.frame(form, data = test_data))
beta_0 <- rep(0, ncol(X))

# Compare fits
for(is_random_offset in c(F, T)){
  offsets <- if(is_random_offset)
    rexp(nrow(test_data), 1) else
      rep(0, nrow(test_data))

  glm_fit <- glm.fit(x = X, y = y, start = beta_0, offset = offsets, family = poisson())


  IWLS_fit <-
    IWLS_poisson(X = X, y = y, beta = rep(0, ncol(X)), offsets = offsets, it_max = 1e2,
               eps = glm.control()$epsilon)

  test_that("glm.fit and IWLS return the same for binomial model",{
    expect_equal(unname(glm_fit$coefficients), c(IWLS_fit))
  })
}
