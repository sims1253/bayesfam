#' Lognormal density distribution in median parametrization.
#'
#' @param x Value space of the distribution, x > 0
#' @param mu Median parameter, mu is already log-transformed, mu unbound
#' @param sigma Sigma shape parameter, sigma >= 0
#' @param log Bool argument, if true, returns the logarithmic density
#'
#' @return Normal distribution density with logit link function
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 10, length.out = 100)
#' plot(x, dlognormal(x, mu = 1, sigma = 0.5), type = "l")
dlognormal <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("lognormal is only defined for x > 0")
  }
  if (isTRUE(sigma <= 0)) {
    stop("lognormal is only defined for sigma > 0")
  }
  # Maybe overkill?
  if (!lenEqual(list_of_vectors = list(x, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("Lognormal_natural argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }
  if (!isLogic_len(log)) {
    stop("the log argument of a density has to be a scalar boolean")
  }

  logpdf <-
    -(log(sigma) + 0.5 * (log(2 * pi))) +
    -log(x) +
    (-(log(x) - mu)^2 / (2 * sigma^2))
  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Lognormal RNG-function in median parametrization.
#'
#' @param n Number of draws
#' @param mu Median parameter, mu unbound, mu already log transformed
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @returns n Lognormal distributed samples
#'
#' @export
#'
#' @examples hist(rlognormal(100, 1, 0.5))
rlognormal <- function(n, mu = 0, sigma = 1) {
  # check the arguments
  if (isTRUE(sigma <= 0)) {
    stop("lognormal is only defined for sigma > 0")
  }
  # Maybe overkill?
  if (!lenEqual(list_of_vectors = list(mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("lognormal argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }
  if (!isNat_len(n)) {
    stop("n has to be a natural scalar number")
  }

  return(
    exp(rnorm(n, mu, sigma))
  )
}
