#' Simulate data from a gamma frailty survival model
#'
#' Generates survival times and censoring indicators under a proportional
#' hazards frailty model with gamma distributed frailty. The baseline hazard can
#' be specified by name or by supplying custom hazard functions.
#'
#' @param baseline Character string identifying the baseline distribution
#'   ("exponential", "weibull", "gompertz", "loglogistic", "lognormal") or a
#'   list with elements `h0`, `H0`, and optionally `H0_inv`, `param_names` and
#'   `positive` indicating which parameters must be positive.
#' @param n Number of observations to simulate.
#' @param params List with elements:
#'   - `baseline`: numeric vector of baseline parameters.
#'   - `beta`: regression coefficients.
#'   - `frailty_var`: variance of the gamma frailty (> 0).
#'   - `censor_rate`: rate for exponential censoring (default 0 for no censoring).
#' @param X Matrix of covariates with `n` rows.
#'
#' @return A list with components `time` and `event`.
#' @examples
#' params <- list(baseline = 0.1, beta = 0.5, frailty_var = 0.2)
#' X <- matrix(rnorm(20), ncol = 1)
#' sim <- simulate_frailty_data("exponential", n = 20, params = params, X = X)
#' init <- list(baseline = 0.1, beta = 0, frailty_var = 0.2)
#' fit <- estimate_frailty_mle(init, sim$time, sim$event, X, "exponential")
#' \dontrun{
#' frailty_simulation_pipeline("exponential", sims = 2, n = c(20, 40),
#'                             true_params = params, plot = FALSE)
#' }
#' @export
simulate_frailty_data <- function(baseline, n, params, X) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  if (!is.numeric(n) || length(n) != 1 || n <= 0) stop("n must be positive")
  if (!is.matrix(X) || nrow(X) != n) stop("X must be a matrix with n rows")
  if (!is.list(params)) stop("params must be a list")
  if (is.null(params$beta)) stop("params$beta is required")
  beta <- params$beta
  p <- ncol(X)
  if (length(beta) != p) stop("length of beta must match ncol(X)")
  frailty_var <- params$frailty_var
  if (!is.numeric(frailty_var) || frailty_var <= 0) {
    stop("frailty_var must be > 0")
  }
  censor_rate <- params$censor_rate %||% 0
  if (!is.numeric(censor_rate) || censor_rate < 0) {
    stop("censor_rate must be >= 0")
  }

  base <- get_baseline_functions(baseline)
  base_params <- params$baseline
  if (length(base$param_names) > 0) {
    if (is.null(base_params) || length(base_params) != length(base$param_names)) {
      stop("baseline parameters missing or incorrect length")
    }
    if (any(base_params[base$positive] <= 0)) {
      stop("baseline parameters violating positivity constraints")
    }
  } else {
    base_params <- numeric()
  }

  Z <- rgamma(n, shape = 1 / frailty_var, scale = frailty_var)
  eta <- as.vector(X %*% beta)
  U <- runif(n)
  time <- base$H0_inv(-log(U) / (Z * exp(eta)), base_params)
  censor <- if (censor_rate > 0) rexp(n, rate = censor_rate) else rep(Inf, n)
  obs_time <- pmin(time, censor)
  event <- as.numeric(time <= censor)
  list(time = obs_time, event = event)
}
