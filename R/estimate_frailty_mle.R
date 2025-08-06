#' Maximum likelihood estimation for gamma frailty model
#'
#' Fits a proportional hazards frailty model with gamma distributed frailty via
#' maximum likelihood. The baseline hazard may be specified by name or through a
#' custom list of hazard functions.
#'
#' @param params_init List with initial parameters:
#'   - `baseline`: numeric vector of baseline parameters.
#'   - `beta`: initial regression coefficients.
#'   - `frailty_var`: initial frailty variance.
#' @param times Observed follow-up times.
#' @param event Event indicators (1 = event, 0 = censored).
#' @param X Matrix of covariates with rows equal to length of `times`.
#' @param baseline Baseline specification (see `simulate_frailty_data`).
#' 
#' @return A list with estimated parameters, standard errors, and log-likelihood.
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
estimate_frailty_mle <- function(params_init, times, event, X, baseline) {
  n <- length(times)
  if (length(event) != n) stop("times and event must have same length")
  if (!is.matrix(X) || nrow(X) != n) stop("X must be a matrix with n rows")

  beta_init <- params_init$beta
  frailty_var_init <- params_init$frailty_var
  base_params_init <- params_init$baseline
  if (is.null(beta_init) || is.null(frailty_var_init)) stop("params_init incomplete")
  p <- ncol(X)
  if (length(beta_init) != p) stop("beta length mismatch with columns of X")

  base <- get_baseline_functions(baseline)
  k <- length(base$param_names)
  if (length(base_params_init) != k) stop("baseline parameter length mismatch")

  init_vec <- c(base_params_init, beta_init, frailty_var_init)
  lower <- c(ifelse(base$positive, 1e-6, -Inf), rep(-Inf, p), 1e-6)

  neg_loglik <- function(par) {
    base_par <- par[1:k]
    beta <- if (p > 0) par[(k + 1):(k + p)] else numeric()
    theta <- par[k + p + 1]
    if (any(base_par[base$positive] <= 0) || theta <= 0) return(Inf)
    eta <- if (p > 0) as.vector(X %*% beta) else rep(0, n)
    h0 <- base$h0(times, base_par)
    H0 <- base$H0(times, base_par)
    if (any(h0 <= 0) || any(H0 < 0)) return(Inf)
    ll_vec <- event * (log(h0) + eta) +
      lgamma(1 / theta + event) - lgamma(1 / theta) -
      (1 / theta + event) * log(1 + theta * H0 * exp(eta)) -
      (1 / theta) * log(theta)
    -sum(ll_vec)
  }

  opt <- optim(init_vec, neg_loglik, method = "L-BFGS-B", lower = lower, hessian = TRUE)
  est <- opt$par
  vc <- try(solve(opt$hessian), silent = TRUE)
  se <- if (inherits(vc, "try-error")) rep(NA, length(est)) else sqrt(diag(vc))
  est_base <- est[1:k]
  est_beta <- if (p > 0) est[(k + 1):(k + p)] else numeric()
  est_theta <- est[k + p + 1]
  se_base <- se[1:k]
  se_beta <- if (p > 0) se[(k + 1):(k + p)] else numeric()
  se_theta <- se[k + p + 1]
  list(
    estimates = list(
      baseline = setNames(est_base, base$param_names),
      beta = est_beta,
      frailty_var = est_theta
    ),
    se = list(
      baseline = setNames(se_base, base$param_names),
      beta = se_beta,
      frailty_var = se_theta
    ),
    logLik = -opt$value
  )
}
