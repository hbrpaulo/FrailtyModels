#' Log-logistic baseline hazard and derivatives
#'
#' Functions for a log-logistic baseline hazard. Parameter vector is
#' `params = c(shape, scale)` with both parameters > 0.
#'
#' @param t Numeric vector of times.
#' @param params Numeric vector `c(shape, scale)`.
#'
#' @return `h0_loglogistic` and `H0_loglogistic` return numeric vectors of hazard
#'   and cumulative hazard values.
#' @examples
#' h0_loglogistic(1:3, c(2, 1))
#' @export
h0_loglogistic <- function(t, params) {
  shape <- params[1]
  scale <- params[2]
  (shape / scale) * (t / scale)^(shape - 1) / (1 + (t / scale)^shape)
}

#' @rdname h0_loglogistic
#' @export
H0_loglogistic <- function(t, params) {
  shape <- params[1]
  scale <- params[2]
  log(1 + (t / scale)^shape)
}

#' Log-likelihood for log-logistic baseline
#'
#' @param params Numeric vector `c(shape, scale)`.
#' @param times Follow-up times.
#' @param event Event indicators (1 = event, 0 = censored).
#'
#' @return Numeric scalar log-likelihood.
#' @export
loglik_loglogistic <- function(params, times, event) {
  shape <- params[1]
  scale <- params[2]
  if (any(params <= 0)) return(-Inf)
  h <- h0_loglogistic(times, params)
  H <- H0_loglogistic(times, params)
  sum(event * log(h) - H)
}

#' Score vector for log-logistic baseline
#'
#' Numerical gradient of the log-likelihood with respect to the parameters.
#'
#' @inheritParams loglik_loglogistic
#'
#' @return Numeric vector of length two.
#' @export
score_loglogistic <- function(params, times, event) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for score_loglogistic")
  }
  numDeriv::grad(loglik_loglogistic, params, times = times, event = event)
}

#' Hessian matrix for log-logistic baseline
#'
#' Numerical Hessian of the log-likelihood with respect to the parameters.
#'
#' @inheritParams loglik_loglogistic
#'
#' @return `2 x 2` numeric matrix.
#' @export
hessian_loglogistic <- function(params, times, event) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for hessian_loglogistic")
  }
  numDeriv::hessian(loglik_loglogistic, params, times = times, event = event)
}

