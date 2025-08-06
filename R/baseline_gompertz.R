#' Gompertz baseline hazard and derivatives
#'
#' Provides hazard, cumulative hazard, log-likelihood, score and Hessian for a
#' Gompertz baseline hazard. Parameter vector is `params = c(lambda, gamma)` with
#' both parameters > 0.
#'
#' @param t Numeric vector of times.
#' @param params Numeric vector `c(lambda, gamma)`.
#'
#' @return `h0_gompertz` and `H0_gompertz` return numeric vectors of hazard and
#'   cumulative hazard values.
#' @examples
#' h0_gompertz(1:3, c(0.5, 0.1))
#' @export
h0_gompertz <- function(t, params) {
  lambda <- params[1]
  gamma <- params[2]
  lambda * exp(gamma * t)
}

#' @rdname h0_gompertz
#' @export
H0_gompertz <- function(t, params) {
  lambda <- params[1]
  gamma <- params[2]
  lambda / gamma * (exp(gamma * t) - 1)
}

#' Log-likelihood for Gompertz baseline
#'
#' @param params Numeric vector `c(lambda, gamma)`.
#' @param times Follow-up times.
#' @param event Event indicators (1 = event, 0 = censored).
#'
#' @return Numeric scalar log-likelihood.
#' @export
loglik_gompertz <- function(params, times, event) {
  lambda <- params[1]
  gamma <- params[2]
  if (any(params <= 0)) return(-Inf)
  h <- h0_gompertz(times, params)
  H <- H0_gompertz(times, params)
  sum(event * log(h) - H)
}

#' Score vector for Gompertz baseline
#'
#' Numerical gradient of the log-likelihood with respect to the parameters.
#'
#' @inheritParams loglik_gompertz
#'
#' @return Numeric vector of length two.
#' @export
score_gompertz <- function(params, times, event) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for score_gompertz")
  }
  numDeriv::grad(loglik_gompertz, params, times = times, event = event)
}

#' Hessian matrix for Gompertz baseline
#'
#' Numerical Hessian of the log-likelihood with respect to the parameters.
#'
#' @inheritParams loglik_gompertz
#'
#' @return `2 x 2` numeric matrix.
#' @export
hessian_gompertz <- function(params, times, event) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for hessian_gompertz")
  }
  numDeriv::hessian(loglik_gompertz, params, times = times, event = event)
}

