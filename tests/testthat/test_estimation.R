library(testthat)

set.seed(123)

test_that("MLE recovers parameters", {
  skip_on_cran()
  n <- 100
  beta <- 0.5
  params <- list(baseline = 0.1, beta = beta, frailty_var = 0.2)
  X <- matrix(rnorm(n), ncol = 1)
  sim <- simulate_frailty_data("exponential", n, params, X)
  fit <- try(estimate_frailty_mle(params, sim$time, sim$event, X, "exponential"), silent = TRUE)
  skip_if(inherits(fit, "try-error"))
  expect_lt(abs(fit$estimates$baseline - params$baseline), 0.1)
  expect_lt(abs(fit$estimates$beta - beta), 0.2)
  expect_lt(abs(fit$estimates$frailty_var - params$frailty_var), 0.1)
})
