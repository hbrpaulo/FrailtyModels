#' Empirical coverage for frailty model parameters
#'
#' Simulates datasets from a gamma frailty model and fits the model to
#' estimate empirical coverage probabilities for nominal confidence
#' intervals.
#'
#' @param baseline Baseline hazard specification (character string or list).
#' @param true_params List with elements `baseline`, `beta`, `frailty_var`,
#'   and optionally `censor_rate`.
#' @param n Number of observations in each simulation.
#' @param sims Number of Monte Carlo replications.
#' @param alpha Significance level for two-sided confidence intervals.
#' @param X_fun Optional function generating an `n x p` covariate matrix;
#'   defaults to standard normal covariates where `p` is the length of
#'   `true_params$beta`.
#' @param save_path Optional file path to save results via `saveRDS`.
#'
#' @return A list with elements
#'   - `coverage`: named numeric vector of empirical coverage probabilities.
#'   - `ci_data`: data frame of confidence interval results across simulations.
#'   If `save_path` is provided, the list is saved to disk.
#' @export
coverage_prob <- function(baseline, true_params, n, sims, alpha = 0.05,
                          X_fun = NULL, save_path = NULL) {
  beta <- true_params$beta
  if (is.null(beta)) beta <- numeric()
  p <- length(beta)
  X_gen <- function(n) {
    if (is.null(X_fun)) {
      if (p > 0) matrix(rnorm(n * p), nrow = n, ncol = p) else matrix(nrow = n, ncol = 0)
    } else {
      X <- X_fun(n)
      if (!is.matrix(X) || nrow(X) != n || ncol(X) != p) {
        stop("X_fun must return an n x p matrix")
      }
      X
    }
  }
  base_names <- names(true_params$baseline)
  if (is.null(base_names)) base_names <- paste0("baseline_", seq_along(true_params$baseline))
  beta_names <- if (p > 0) paste0("beta_", seq_len(p)) else character(0)
  param_names <- c(base_names, beta_names, "frailty_var")
  n_par <- length(param_names)
  coverage_mat <- matrix(NA, nrow = sims, ncol = n_par)
  colnames(coverage_mat) <- param_names
  ci_records <- vector("list", sims)
  true_vals <- c(true_params$baseline, beta, true_params$frailty_var)
  names(true_vals) <- param_names
  for (i in seq_len(sims)) {
    X <- X_gen(n)
    sim <- simulate_frailty_data(baseline, n, params = true_params, X = X)
    params_init <- true_params
    params_init$censor_rate <- NULL
    fit <- try(estimate_frailty_mle(params_init, sim$time, sim$event, X, baseline),
               silent = TRUE)
    if (inherits(fit, "try-error")) next
    est_base <- fit$estimates$baseline
    est_beta <- fit$estimates$beta
    est_theta <- fit$estimates$frailty_var
    se_base <- fit$se$baseline
    se_beta <- fit$se$beta
    se_theta <- fit$se$frailty_var
    est <- c(est_base, est_beta, est_theta)
    se <- c(se_base, se_beta, se_theta)
    names(est) <- param_names
    z <- qnorm(1 - alpha / 2)
    lower <- est - z * se
    upper <- est + z * se
    coverage_mat[i, ] <- (true_vals >= lower) & (true_vals <= upper)
    ci_records[[i]] <- data.frame(
      Parameter = param_names,
      Lower = lower,
      Estimate = est,
      Upper = upper,
      True = true_vals,
      Sim = i
    )
  }
  keep <- complete.cases(coverage_mat)
  coverage <- colMeans(coverage_mat[keep, , drop = FALSE])
  ci_all <- do.call(rbind, ci_records[keep])
  res <- list(coverage = coverage, ci_data = ci_all)
  if (!is.null(save_path)) saveRDS(res, save_path)
  res
}
