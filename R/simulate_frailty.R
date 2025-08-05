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
      if (length(baseline) == 1 && is.character(baseline[[1]])) {
        return(resolve_baseline(baseline[[1]]))
      }
      required <- c("h0", "H0")
      if (!all(required %in% names(baseline))) {
        stop("Custom baseline must provide h0 and H0 functions")
      }
      if (!is.function(baseline$h0) || !is.function(baseline$H0)) {
        stop("h0 and H0 must be functions")
      }
      b <- baseline
      b$param_names <- b$param_names %||% character(0)
      b$positive <- b$positive %||% rep(TRUE, length(b$param_names))
      if (is.null(b$H0_inv)) {
        H0_fun <- b$H0
        b$H0_inv <- function(y, p) {
          sapply(y, function(yy) {
            uniroot(function(t) H0_fun(t, p) - yy,
                    lower = 1e-10, upper = 1e6)$root
          })
        }
      }
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
