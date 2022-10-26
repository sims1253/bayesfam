#' Mean parameterization the beta pdf.
#'
#' @param x x-value, x e (0, 1)
#' @param mu Mean parameter, mu e (0, 1)
#' @param phi Precision parameter, phi > 0
#' @param log Optional argument. If TRUE, returns log(pdf). Normally False.
#'
#' @details The Beta Distribution has Density
#' \deqn{f(y | \mu, \phi) = \frac{\Gamma(\phi) x^{\mu\phi - 1} (1 - x)^{(1-\mu)\phi}}{\Gamma(\mu\phi)\Gamma((1 - \mu)\phi)} }
#' @details With parameterisation of the usual Beta-Distribution's shape parameters a and b as:
#' \deqn{a := \mu\phi, b := (1 - \mu)\phi}
#'
#' @return PDF of custom Beta Distribution
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, dbeta_mean(x, mu = 0.5, phi = 1), type = "l")
dbeta_mean <- function(x, mu, phi, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("The value x has to be in (0, 1).")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0, 1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  lpdf <- (log(gamma(phi)) -
    log(gamma(mu * phi)) -
    log(gamma((1 - mu) * phi))) +
    log(x) * (mu * phi - 1) + log1p(-x) * ((1 - mu) * phi - 1)
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Mean parameterization of the beta quantile function
#'
#' @param p Probabilities, for which to calculate the quantiles
#' @param mu Mean, mu e (0, 1).
#' @param phi Precision parameter, phi > 0.
#'
#' @return Quantiles of the beta distribution
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' plot(x, qbeta_mean(x, mu = 0.5, phi = 2), type = "l")
qbeta_mean <- function(p, mu, phi) {
  if (isTRUE(any(p < 0 | p > 1))) {
    stop("p has to be in an interval of [0, 1]")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0, 1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  qbeta(p, shape1 = mu * phi, shape2 = (1 - mu) * phi)
}

#' Beta distribution RNG for mean parameterization
#'
#' @param n Number of draws.
#' @param mu Mean
#' @param phi Precision
#'
#' @return n samples Beta distributed.
#' @export
#' @examples hist(rbeta_mean(1000, mu = 0.5, phi = 1))
rbeta_mean <- function(n, mu, phi) {
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0,1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  return(rbeta(n, mu * phi, (1 - mu) * phi))
}
