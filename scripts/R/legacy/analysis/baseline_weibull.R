#' Weibull baseline hazard function
#'
#' Computes the baseline hazard for a Weibull distribution.
#'
#' @param t Numeric vector of times.
#' @param shape Weibull shape parameter (>0).
#' @param scale Weibull scale parameter (>0).
#'
#' @return Numeric vector of baseline hazards evaluated at `t`.
#' @export
h0_weibull <- function(t, shape, scale) {
  (shape / scale) * (t / scale)^(shape - 1)
}

#' Weibull cumulative baseline hazard function
#'
#' Computes the cumulative baseline hazard for a Weibull distribution.
#'
#' @param t Numeric vector of times.
#' @param shape Weibull shape parameter (>0).
#' @param scale Weibull scale parameter (>0).
#'
#' @return Numeric vector of cumulative baseline hazards evaluated at `t`.
#' @export
H0_weibull <- function(t, shape, scale) {
  (t / scale)^shape
}
