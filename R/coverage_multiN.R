#' Coverage summary across sample sizes
#'
#' Runs `coverage_prob()` for a vector of sample sizes and summarises the
#' resulting coverage probabilities and average confidence interval widths.
#'
#' @param baseline Baseline hazard specification (character string or list).
#' @param true_params List with true parameter values.
#' @param n_vec Integer vector of sample sizes to evaluate.
#' @param sims Number of simulations per sample size.
#' @param alpha Significance level for two-sided confidence intervals.
#' @param X_fun Optional function generating covariate matrices.
#' @param save_path Optional file path to save aggregated results via `saveRDS`.
#'
#' @return A list with elements
#'   - `coverage`: data frame with columns `N`, `Parameter`, and `Coverage`.
#'   - `ci_width`: data frame with columns `N`, `Parameter`, and `MeanWidth`.
#'   If `save_path` is supplied the list is saved to disk.
#' @export
coverage_multiN <- function(baseline, true_params, n_vec, sims, alpha = 0.05,
                            X_fun = NULL, save_path = NULL) {
  cov_records <- list()
  width_records <- list()
  for (n in n_vec) {
    res <- coverage_prob(baseline, true_params, n, sims, alpha = alpha,
                         X_fun = X_fun)
    cov_vec <- res$coverage
    cov_records[[length(cov_records) + 1]] <- data.frame(
      N = n,
      Parameter = names(cov_vec),
      Coverage = as.numeric(cov_vec)
    )
    ci_df <- res$ci_data
    if (is.null(ci_df) || nrow(ci_df) == 0) {
      width_summary <- data.frame(Parameter = names(cov_vec), Width = NA_real_)
    } else {
      ci_df$Width <- ci_df$Upper - ci_df$Lower
      width_summary <- aggregate(Width ~ Parameter, data = ci_df, FUN = mean)
    }
    width_records[[length(width_records) + 1]] <- data.frame(
      N = n,
      Parameter = width_summary$Parameter,
      MeanWidth = width_summary$Width
    )
  }
  coverage_df <- do.call(rbind, cov_records)
  width_df <- do.call(rbind, width_records)
  result <- list(coverage = coverage_df, ci_width = width_df)
  if (!is.null(save_path)) saveRDS(result, save_path)
  result
}
