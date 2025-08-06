#' Frailty simulation pipeline for Weibull baseline
#'
#' Runs a complete simulation pipeline: simulate one data set,
#' estimate parameters, perform coverage studies for one or more
#' sample sizes, and produce a funnel plot of confidence interval widths.
#'
#' @param shape True Weibull shape parameter.
#' @param scale True Weibull scale parameter.
#' @param theta True frailty variance.
#' @param nsim Number of simulations for coverage studies.
#' @param sample_sizes Integer vector of sample sizes to study. The first
#'   element is used for the example simulation and estimation output.
#' @param conf_level Confidence level for intervals.
#' @param censor_rate Rate for exponential censoring used in simulation.
#'
#' @return A list with elements `simulation` (last simulated data set),
#'   `fit` (MLE results), `coverage` (data frame of coverage results), and
#'   `funnel_plot` (ggplot object).
#' @export
frailty_simulation_pipeline_weibull <- function(shape, scale, theta,
                                               nsim = 100,
                                               sample_sizes,
                                               conf_level = 0.95,
                                               censor_rate = 0) {
  sim <- simulate_frailty(sample_sizes[1], shape, scale, theta, censor_rate)
  fit <- estimate_frailty_mle(sim)
  coverage <- coverage_study_multi(sample_sizes, shape, scale, theta,
                                   nsim, conf_level, censor_rate)
  plot <- funnel_plot(coverage)

  list(simulation = sim, fit = fit, coverage = coverage, funnel_plot = plot)
}
