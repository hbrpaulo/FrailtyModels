#' Weibull baseline hazard and derivatives
#'
#' Provides hazard, cumulative hazard, log-likelihood, score, and Hessian for
#' a Weibull baseline hazard. Parameter vector is `params = c(shape, scale)`
#' with both parameters > 0.
#'
#' @param t Numeric vector of times.
#' @param params Numeric vector `c(shape, scale)`.
#'
#' @return `h0_weibull` and `H0_weibull` return numeric vectors with hazard and
#'   cumulative hazard values.
#' @examples
#' h0_weibull(1:3, c(2, 1))
#' @export
h0_weibull <- function(t, params) {
  shape <- params[1]
  scale <- params[2]
  shape / scale * (t / scale)^(shape - 1)
}

#' @rdname h0_weibull
#' @export
H0_weibull <- function(t, params) {
  shape <- params[1]
  scale <- params[2]
  (t / scale)^shape
}

#' Log-likelihood for Weibull baseline
#'
#' @param params Numeric vector `c(shape, scale)`.
#' @param times Follow-up times.
#' @param event Event indicators (1 = event, 0 = censored).
#'
#' @return Numeric scalar log-likelihood.
#' @export
loglik_weibull <- function(params, times, event) {
  shape <- params[1]
  scale <- params[2]
  if (any(params <= 0)) return(-Inf)
  h <- h0_weibull(times, params)
  H <- H0_weibull(times, params)
  sum(event * log(h) - H)
}

#' Score vector for Weibull baseline
#'
#' @inheritParams loglik_weibull
#'
#' @return Numeric vector of length two containing the score for `shape` and
#'   `scale`.
#' @export
score_weibull <- function(params, times, event) {
  shape <- params[1]
  scale <- params[2]
  z <- times / scale
  log_z <- log(z)
  H <- z^shape
  s1 <- sum(event * (1 / shape + log_z) - H * log_z)
  s2 <- sum((shape / scale) * (H - event))
  c(s1, s2)
}

#' Hessian matrix for Weibull baseline
#'
#' @inheritParams loglik_weibull
#'
#' @return `2 x 2` numeric matrix with second derivatives of the
#'   log-likelihood.
#' @export
hessian_weibull <- function(params, times, event) {
  shape <- params[1]
  scale <- params[2]
  z <- times / scale
  log_z <- log(z)
  H <- z^shape
  h11 <- sum(-event / shape^2 - H * log_z^2)
  h12 <- sum(-event / scale + H * (shape * log_z + 1) / scale)
  h22 <- sum((shape / scale^2) * (event - (1 + shape) * H))
  matrix(c(h11, h12, h12, h22), nrow = 2, byrow = TRUE)
}
