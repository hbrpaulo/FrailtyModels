#' Log-normal baseline hazard and derivatives
#'
#' Hazard-related functions for a log-normal baseline distribution. Parameter
#' vector is `params = c(meanlog, sdlog)` where `sdlog > 0`.
#'
#' @param t Numeric vector of times.
#' @param params Numeric vector `c(meanlog, sdlog)`.
#'
#' @return `h0_lognormal` and `H0_lognormal` return numeric vectors of hazard and
#'   cumulative hazard values.
#' @examples
#' h0_lognormal(1:3, c(0, 1))
#' @export
h0_lognormal <- function(t, params) {
  meanlog <- params[1]
  sdlog <- params[2]
  dlnorm(t, meanlog = meanlog, sdlog = sdlog) /
    (1 - plnorm(t, meanlog = meanlog, sdlog = sdlog))
}

#' @rdname h0_lognormal
#' @export
H0_lognormal <- function(t, params) {
  meanlog <- params[1]
  sdlog <- params[2]
  -log(1 - plnorm(t, meanlog = meanlog, sdlog = sdlog))
}

#' Log-likelihood for log-normal baseline
#'
#' @param params Numeric vector `c(meanlog, sdlog)`.
#' @param times Follow-up times.
#' @param event Event indicators (1 = event, 0 = censored).
#'
#' @return Numeric scalar log-likelihood.
#' @export
loglik_lognormal <- function(params, times, event) {
  meanlog <- params[1]
  sdlog <- params[2]
  if (sdlog <= 0) return(-Inf)
  h <- h0_lognormal(times, params)
  H <- H0_lognormal(times, params)
  sum(event * log(h) - H)
}

#' Score vector for log-normal baseline
#'
#' Numerical gradient of the log-likelihood with respect to the parameters.
#'
#' @inheritParams loglik_lognormal
#'
#' @return Numeric vector of length two.
#' @export
score_lognormal <- function(params, times, event) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for score_lognormal")
  }
  numDeriv::grad(loglik_lognormal, params, times = times, event = event)
}

#' Hessian matrix for log-normal baseline
#'
#' Numerical Hessian of the log-likelihood with respect to the parameters.
#'
#' @inheritParams loglik_lognormal
#'
#' @return `2 x 2` numeric matrix.
#' @export
hessian_lognormal <- function(params, times, event) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required for hessian_lognormal")
  }
  numDeriv::hessian(loglik_lognormal, params, times = times, event = event)
}
