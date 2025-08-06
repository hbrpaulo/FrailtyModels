#' Retrieve baseline hazard functions
#'
#' Resolves baseline hazard functions either from a character string naming a
#' supported distribution or from a custom list of functions.
#'
#' @param baseline Either a character string ("exponential", "weibull",
#'   "gompertz", "loglogistic", "lognormal") or a list with elements `h0` and
#'   `H0`, and optionally `H0_inv`, `param_names`, and `positive` indicating which
#'   parameters must be positive.
#'
#' @return A list with components `h0`, `H0`, `H0_inv`, `param_names` and
#'   `positive`.
#' @examples
#' get_baseline_functions("exponential")
#' @export
get_baseline_functions <- function(baseline) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  if (is.character(baseline)) {
    b <- switch(tolower(baseline),
      exponential = list(
        h0 = h0_exponential,
        H0 = H0_exponential,
        H0_inv = function(y, p) y / p[1],
        param_names = "rate",
        positive = TRUE
      ),
      weibull = list(
        h0 = h0_weibull,
        H0 = H0_weibull,
        H0_inv = function(y, p) p[2] * y^(1 / p[1]),
        param_names = c("shape", "scale"),
        positive = c(TRUE, TRUE)
      ),
      gompertz = list(
        h0 = h0_gompertz,
        H0 = H0_gompertz,
        H0_inv = function(y, p) log(1 + p[2] * y / p[1]) / p[2],
        param_names = c("lambda", "gamma"),
        positive = c(TRUE, TRUE)
      ),
      loglogistic = list(
        h0 = h0_loglogistic,
        H0 = H0_loglogistic,
        H0_inv = function(y, p) p[2] * (exp(y) - 1)^(1 / p[1]),
        param_names = c("shape", "scale"),
        positive = c(TRUE, TRUE)
      ),
      lognormal = list(
        h0 = h0_lognormal,
        H0 = H0_lognormal,
        H0_inv = function(y, p) qlnorm(1 - exp(-y), meanlog = p[1], sdlog = p[2]),
        param_names = c("meanlog", "sdlog"),
        positive = c(FALSE, TRUE)
      ),
      stop("Unsupported baseline distribution")
    )
  } else if (is.list(baseline)) {
    if (length(baseline) == 1 && is.character(baseline[[1]])) {
      return(get_baseline_functions(baseline[[1]]))
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
