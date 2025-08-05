#' Piecewise exponential baseline hazard and derivatives
#'
#' Functions for a piecewise constant baseline hazard. Parameters are the vector
#' of rates for each interval defined by `breaks`.
#'
#' @param t Numeric vector of times.
#' @param breaks Increasing numeric vector of interval boundaries of length
#'   `length(rates) + 1`.
#' @param rates Numeric vector of interval rates.
#'
#' @return `h0_piecewise_exp` and `H0_piecewise_exp` return numeric vectors of
#'   hazard and cumulative hazard values.
#' @examples
#' h0_piecewise_exp(c(1, 3, 6), c(0, 2, 5, 10), c(0.5, 1, 0.2))
#' @export
h0_piecewise_exp <- function(t, breaks, rates) {
  rates[findInterval(t, breaks, rightmost.closed = TRUE)]
}

#' @rdname h0_piecewise_exp
#' @export
H0_piecewise_exp <- function(t, breaks, rates) {
  sapply(t, function(x) {
    idx <- findInterval(x, breaks, rightmost.closed = TRUE)
    if (idx == 0) return(0)
    cum <- sum(rates[seq_len(idx - 1)] * diff(breaks[seq_len(idx)]))
    cum + rates[idx] * (x - breaks[idx])
  })
}

#' Log-likelihood for piecewise exponential baseline
#'
#' @param rates Numeric vector of interval rates.
#' @param times Follow-up times.
#' @param event Event indicators (1 = event, 0 = censored).
#' @param breaks Increasing numeric vector of interval boundaries.
#'
#' @return Numeric scalar log-likelihood.
#' @export
loglik_piecewise_exp <- function(rates, times, event, breaks) {
  if (any(rates <= 0)) return(-Inf)
  interval <- findInterval(times, breaks, rightmost.closed = TRUE)
  exposure <- sapply(seq_along(rates), function(k) {
    pmax(pmin(times, breaks[k + 1]) - breaks[k], 0)
  })
  event_counts <- tabulate(interval[event == 1], nbins = length(rates))
  E_k <- colSums(exposure)
  sum(event_counts * log(rates) - rates * E_k)
}

#' Score vector for piecewise exponential baseline
#'
#' @inheritParams loglik_piecewise_exp
#'
#' @return Numeric vector of scores for each rate.
#' @export
score_piecewise_exp <- function(rates, times, event, breaks) {
  interval <- findInterval(times, breaks, rightmost.closed = TRUE)
  exposure <- sapply(seq_along(rates), function(k) {
    pmax(pmin(times, breaks[k + 1]) - breaks[k], 0)
  })
  event_counts <- tabulate(interval[event == 1], nbins = length(rates))
  E_k <- colSums(exposure)
  event_counts / rates - E_k
}

#' Hessian matrix for piecewise exponential baseline
#'
#' @inheritParams loglik_piecewise_exp
#'
#' @return Diagonal matrix with second derivatives of the log-likelihood.
#' @export
hessian_piecewise_exp <- function(rates, times, event, breaks) {
  interval <- findInterval(times, breaks, rightmost.closed = TRUE)
  event_counts <- tabulate(interval[event == 1], nbins = length(rates))
  diag(-event_counts / rates^2)
}
