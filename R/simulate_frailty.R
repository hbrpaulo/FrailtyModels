#' Simulate data from a gamma frailty survival model
#'
#' Generates survival times and censoring indicators under a proportional
#' hazards frailty model with gamma distributed frailty. The baseline hazard can
#' be specified by name or by supplying custom hazard functions.
#'
#' @param baseline Character string identifying the baseline distribution
#'   ("exponential", "weibull", "gompertz", "loglogistic", "lognormal") or a
#'   list with elements `h0`, `H0`, `H0_inv`, and optionally `param_names` and
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
#' @export
simulate_frailty_data <- function(baseline, n, params, X) {
  resolve_baseline <- function(baseline) {
    if (is.character(baseline)) {
      b <- switch(tolower(baseline),
        exponential = list(
          param_names = "rate",
          positive = TRUE,
          h0 = function(t, p) rep(p[1], length(t)),
          H0 = function(t, p) p[1] * t,
          H0_inv = function(y, p) y / p[1]
        ),
        weibull = list(
          param_names = c("shape", "scale"),
          positive = c(TRUE, TRUE),
          h0 = function(t, p) p[1] / p[2] * (t / p[2])^(p[1] - 1),
          H0 = function(t, p) (t / p[2])^p[1],
          H0_inv = function(y, p) p[2] * y^(1 / p[1])
        ),
        gompertz = list(
          param_names = c("lambda", "gamma"),
          positive = c(TRUE, TRUE),
          h0 = function(t, p) p[1] * exp(p[2] * t),
          H0 = function(t, p) p[1] / p[2] * (exp(p[2] * t) - 1),
          H0_inv = function(y, p) log(1 + p[2] * y / p[1]) / p[2]
        ),
        loglogistic = list(
          param_names = c("shape", "scale"),
          positive = c(TRUE, TRUE),
          h0 = function(t, p) (p[1] / p[2]) * (t / p[2])^(p[1] - 1) /
            (1 + (t / p[2])^p[1]),
          H0 = function(t, p) log(1 + (t / p[2])^p[1]),
          H0_inv = function(y, p) p[2] * (exp(y) - 1)^(1 / p[1])
        ),
        lognormal = list(
          param_names = c("meanlog", "sdlog"),
          positive = c(FALSE, TRUE),
          h0 = function(t, p) dlnorm(t, meanlog = p[1], sdlog = p[2]) /
            (1 - plnorm(t, meanlog = p[1], sdlog = p[2])),
          H0 = function(t, p) -log(1 - plnorm(t, meanlog = p[1], sdlog = p[2])),
          H0_inv = function(y, p) qlnorm(1 - exp(-y), meanlog = p[1], sdlog = p[2])
        ),
        stop("Unsupported baseline distribution")
      )
    } else if (is.list(baseline)) {
      required <- c("h0", "H0", "H0_inv")
      if (!all(required %in% names(baseline))) {
        stop("Custom baseline must provide h0, H0, and H0_inv functions")
      }
      b <- baseline
      b$param_names <- b$param_names %||% character(0)
      b$positive <- b$positive %||% rep(TRUE, length(b$param_names))
    } else {
      stop("baseline must be a character string or list")
    }
    b
  }
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

  base <- resolve_baseline(baseline)
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
