#' Frailty model simulation pipeline
#'
#' Runs repeated simulations for a gamma frailty model, fits the model by
#' maximum likelihood, assesses confidence interval coverage and optionally
#' draws a funnel plot of interval width against sample size.
#'
#' @param baseline Either a character string identifying the baseline
#'   distribution ("exponential", "weibull", "gompertz", "loglogistic" or
#'   "lognormal") or a list with elements `h0`, `H0` and optionally `H0_inv`.
#' @param sims Number of simulation replicates.
#' @param n Integer vector of sample sizes.
#' @param true_params List with components `baseline` (numeric vector of true
#'   baseline parameters), `beta` (regression coefficients), `frailty_var`
#'   (gamma frailty variance) and optionally `censor_rate` (rate for exponential
#'   censoring).
#' @param plot Logical; if `TRUE`, produce a funnel plot of confidence interval
#'   width versus sample size.
#'
#' @return A list with elements `coverage` (data frame of coverage
#'   probabilities), `ci_width` (data frame of confidence interval widths) and
#'   `plot` (ggplot object when `plot = TRUE`).
#'
#' @examples
#' params <- list(baseline = 0.1, beta = 0.5, frailty_var = 0.2)
#' X <- matrix(rnorm(20), ncol = 1)
#' sim <- simulate_frailty_data("exponential", n = 20, params = params, X = X)
#' fit <- estimate_frailty_mle(list(baseline = 0.1, beta = 0,
#'                                  frailty_var = 0.2),
#'                             sim$time, sim$event, X, "exponential")
#' \dontrun{
#' frailty_simulation_pipeline("exponential", sims = 2, n = c(20, 40),
#'                             true_params = params, plot = FALSE)
#' }
#' @export
frailty_simulation_pipeline <- function(baseline, sims, n, true_params, plot = TRUE) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  base <- get_baseline_functions(baseline)
  n_vec <- n
  p <- length(true_params$beta)
  param_names <- c(base$param_names,
                   if (p > 0) paste0("beta", seq_len(p)) else character(0),
                   "frailty_var")
  true_vec <- c(true_params$baseline, true_params$beta, true_params$frailty_var)

  coverage_list <- list()
  ci_records <- list()
  for (n_i in n_vec) {
    cover_mat <- matrix(NA, nrow = sims, ncol = length(true_vec))
    width_mat <- matrix(NA, nrow = sims, ncol = length(true_vec))
    for (s in seq_len(sims)) {
      X <- if (p > 0) matrix(rnorm(n_i * p), nrow = n_i, ncol = p) else matrix(nrow = n_i, ncol = 0)
      sim <- simulate_frailty_data(
        baseline = base,
        n = n_i,
        params = list(
          baseline = true_params$baseline,
          beta = true_params$beta,
          frailty_var = true_params$frailty_var,
          censor_rate = true_params$censor_rate %||% 0
        ),
        X = X
      )
      init <- list(
        baseline = true_params$baseline,
        beta = rep(0, p),
        frailty_var = true_params$frailty_var
      )
      fit <- try(estimate_frailty_mle(init, sim$time, sim$event, X, base), silent = TRUE)
      if (inherits(fit, "try-error")) next
      est <- c(fit$estimates$baseline, fit$estimates$beta, fit$estimates$frailty_var)
      se <- c(fit$se$baseline, fit$se$beta, fit$se$frailty_var)
      if (any(is.na(se))) next
      z <- qnorm(0.975)
      lower <- est - z * se
      upper <- est + z * se
      cover_mat[s, ] <- (true_vec >= lower) & (true_vec <= upper)
      width_mat[s, ] <- upper - lower
    }
    coverage <- colMeans(cover_mat, na.rm = TRUE)
    coverage_list[[as.character(n_i)]] <- data.frame(
      n = n_i,
      parameter = param_names,
      coverage = coverage
    )
    if (sum(!is.na(width_mat)) > 0) {
      ci_records[[length(ci_records) + 1]] <- data.frame(
        width = as.vector(width_mat),
        parameter = rep(param_names, each = nrow(width_mat)),
        n = n_i
      )
    }
  }
  coverage_df <- do.call(rbind, coverage_list)
  ci_width_df <- if (length(ci_records) > 0) do.call(rbind, ci_records) else data.frame()
  plot_obj <- NULL
  if (plot && nrow(ci_width_df) > 0) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' is required for plotting")
    } else {
      plot_obj <- ggplot2::ggplot(ci_width_df, ggplot2::aes(x = n, y = width, color = parameter)) +
        ggplot2::geom_point(alpha = 0.4) +
        ggplot2::geom_smooth(se = FALSE) +
        ggplot2::labs(x = "Sample size", y = "CI width", title = "Funnel plot of CI width")
      print(plot_obj)
    }
  }
  list(coverage = coverage_df, ci_width = ci_width_df, plot = plot_obj)
}

