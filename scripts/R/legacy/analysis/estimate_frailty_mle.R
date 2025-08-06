#' Maximum likelihood estimation for Weibull baseline with Gamma frailty
#'
#' Estimates the Weibull shape and scale parameters and the Gamma frailty
#' variance using maximum likelihood.
#'
#' @param data A data frame with columns `time` and `event` as produced by
#'   [simulate_frailty()].
#' @param start Optional named vector of starting values on the natural
#'   scale for `shape`, `scale`, and `theta`. Defaults are 1.
#'
#' @return A list containing parameter estimates (`coefficients`),
#'   standard errors (`se`), the observed information matrix (`vcov`), and
#'   the maximized log-likelihood (`loglik`).
#' @export
estimate_frailty_mle <- function(data, start = c(shape = 1, scale = 1, theta = 0.5)) {
  neg_loglik <- function(par) {
    shape <- exp(par[1])
    scale <- exp(par[2])
    theta <- exp(par[3])

    H <- (data$time / scale)^shape
    h <- (shape / scale) * (data$time / scale)^(shape - 1)

    loglik <- data$event * (log(h) - (1/theta + 1) * log1p(theta * H)) +
      (1 - data$event) * (-(1/theta) * log1p(theta * H))
    -sum(loglik)
  }

  fit <- stats::optim(log(start), neg_loglik, hessian = TRUE)
  est <- exp(fit$par)
  vcov <- solve(fit$hessian)
  se <- sqrt(diag(vcov)) * est
  loglik <- -fit$value

  list(
    coefficients = stats::setNames(est, c("shape", "scale", "theta")),
    se = stats::setNames(se, c("shape", "scale", "theta")),
    vcov = vcov,
    loglik = loglik
  )
}
