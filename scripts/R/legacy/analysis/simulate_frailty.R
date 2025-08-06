#' Simulate survival times with Gamma frailty and Weibull baseline
#'
#' Generates survival times under a Weibull baseline hazard with
#' individual-level Gamma frailty and optional exponential censoring.
#'
#' @param n Number of individuals to simulate.
#' @param shape Weibull shape parameter.
#' @param scale Weibull scale parameter.
#' @param theta Variance of the Gamma frailty distribution. Frailties are
#'   distributed as `Gamma(1/theta, 1/theta)` giving mean 1 and variance theta.
#' @param censor_rate Rate parameter of the exponential censoring
#'   distribution. Use 0 for no censoring.
#'
#' @return A data frame with columns `time`, `event` (1 = observed event,
#'   0 = censored) and `frailty` (individual frailty values).
#' @export
simulate_frailty <- function(n, shape, scale, theta, censor_rate = 0) {
  w <- rgamma(n, shape = 1/theta, rate = 1/theta)
  u <- runif(n)
  t <- scale * (-log(u) / w)^(1 / shape)

  if (censor_rate > 0) {
    c <- rexp(n, rate = censor_rate)
    event <- as.integer(t <= c)
    time <- pmin(t, c)
  } else {
    event <- rep(1L, n)
    time <- t
  }

  data.frame(time = time, event = event, frailty = w)
}
