#' Exponential baseline hazard and derivatives
#'
#' Functions related to a parametric survival model with an exponential
#' baseline hazard. Parameter vector is `params = c(rate)` with rate > 0.
#'
#' @param t Numeric vector of times.
#' @param params Numeric vector of length 1 giving the rate.
#'
#' @return `h0_exponential` and `H0_exponential` return numeric vectors of
#'   hazard and cumulative hazard values, respectively.
#' @examples
#' h0_exponential(1:3, 0.5)
#' @export
h0_exponential <- function(t, params) {
  rate <- params[1]
  rep(rate, length(t))
}

#' @rdname h0_exponential
#' @export
H0_exponential <- function(t, params) {
  rate <- params[1]
  rate * t
}

#' Log-likelihood for exponential baseline
#'
#' Computes the log-likelihood for independent survival times with an
#' exponential baseline hazard.
#'
#' @param params Numeric vector `c(rate)` with rate > 0.
#' @param times Follow-up times.
#' @param event Event indicators (1 = event, 0 = censored).
#'
#' @return Numeric scalar log-likelihood.
#' @export
loglik_exponential <- function(params, times, event) {
  rate <- params[1]
  if (rate <= 0) return(-Inf)
  sum(event * log(rate) - rate * times)
}

#' Score vector for exponential baseline
#'
#' @inheritParams loglik_exponential
#'
#' @return Numeric vector of length 1 with the score.
#' @export
score_exponential <- function(params, times, event) {
  rate <- params[1]
  c(sum(event / rate - times))
}

#' Hessian matrix for exponential baseline
#'
#' @inheritParams loglik_exponential
#'
#' @return `1 x 1` numeric matrix with the second derivative of the
#'   log-likelihood.
#' @export
hessian_exponential <- function(params, times, event) {
  rate <- params[1]
  matrix(-sum(event) / rate^2, nrow = 1, ncol = 1)
}
