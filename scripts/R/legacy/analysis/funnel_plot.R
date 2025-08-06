#' Funnel plot of confidence interval widths by sample size
#'
#' Produces a plot of mean confidence interval widths versus sample size
#' from the results of [coverage_study_multi()].
#'
#' @param coverage_results Data frame returned by
#'   [coverage_study_multi()].
#'
#' @return A `ggplot2` object visualizing mean CI width against sample
#'   size, faceted by parameter.
#' @export
funnel_plot <- function(coverage_results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for funnel_plot().")
  }

  ggplot2::ggplot(coverage_results, ggplot2::aes(x = n, y = mean_ci_width)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~parameter, scales = "free_y") +
    ggplot2::labs(x = "Sample size", y = "Mean CI width",
                  title = "Funnel plot of CI width by sample size") +
    ggplot2::theme_minimal()
}
