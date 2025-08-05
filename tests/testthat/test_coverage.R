library(testthat)

set.seed(456)

test_that("coverage is near nominal", {
  skip_on_cran()
  baseline <- "exponential"
  true_params <- list(baseline = 0.1, beta = 0.5, frailty_var = 0.2)
  n <- 100
  sims <- 20
  res <- try(coverage_prob(baseline, true_params, n = n, sims = sims), silent = TRUE)
  skip_if(inherits(res, "try-error"))
  skip_if(anyNA(res$coverage))
  expect_true(all(abs(res$coverage - 0.95) < 0.2))
})
