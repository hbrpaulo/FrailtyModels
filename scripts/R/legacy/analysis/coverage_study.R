#' Coverage study for Weibull Gamma frailty model
#'
#' Conducts a coverage study for parameter estimates of the Weibull
#' baseline with Gamma frailty for a fixed sample size.
#'
#' @param n Sample size.
#' @param shape True Weibull shape parameter.
#' @param scale True Weibull scale parameter.
#' @param theta True frailty variance.
#' @param nsim Number of simulation replicates.
#' @param conf_level Confidence level for intervals.
#' @param censor_rate Rate of exponential censoring. Use 0 for none.
#'
#' @return A list with elements `coverage` and `mean_ci_width` (named
#'   vectors for parameters `shape`, `scale`, `theta`).
#' @export
coverage_study <- function(n, shape, scale, theta, nsim = 100,
                           conf_level = 0.95, censor_rate = 0) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  cover_counts <- c(shape = 0, scale = 0, theta = 0)
  width_sums <- c(shape = 0, scale = 0, theta = 0)
  truth <- c(shape = shape, scale = scale, theta = theta)

  for (i in seq_len(nsim)) {
    dat <- simulate_frailty(n, shape, scale, theta, censor_rate)
    fit <- estimate_frailty_mle(dat)
    est <- fit$coefficients
    se <- fit$se
    lower <- est - z * se
    upper <- est + z * se
    cover_counts <- cover_counts + (lower <= truth & upper >= truth)
    width_sums <- width_sums + (upper - lower)
  }

  list(
    coverage = cover_counts / nsim,
    mean_ci_width = width_sums / nsim
  )
}

#' Coverage study over multiple sample sizes
#'
#' Repeats [coverage_study()] over a vector of sample sizes.
#'
#' @param sample_sizes Numeric vector of sample sizes.
#' @inheritParams coverage_study
#'
#' @return A data frame with columns `n`, `parameter`, `coverage`, and
#'   `mean_ci_width`.
#' @export
coverage_study_multi <- function(sample_sizes, shape, scale, theta,
                                 nsim = 100, conf_level = 0.95,
                                 censor_rate = 0) {
  res <- lapply(sample_sizes, function(n) {
    out <- coverage_study(n, shape, scale, theta, nsim, conf_level,
                          censor_rate)
    data.frame(
      n = n,
      parameter = names(out$coverage),
      coverage = as.numeric(out$coverage),
      mean_ci_width = as.numeric(out$mean_ci_width)
    )
  })
  do.call(rbind, res)
}
