#' Mean parameterization of the gamma pdf.
#'
#' @param x Value space, x > 0.
#' @param mu Mean parameter of the density, mu > 0.
#' @param a Shape parameter, a > 0.
#' @param log optional argument. If true, returns logarathmic probability. Default = FALSE
#'
#' @details Define rate constante rho as:
#' \deqn{\rho(\alpha, \mu) = \frac{\alpha}{\mu}}
#' @details The Frechet distribution density is defined as
#' \deqn{f(y) = \frac{y^{\alpha - 1} exp(-\frac{y}{\rho})} {\rho^\alpha \Gamma(\alpha)} }
#'
#' @return f(x | mu, k)
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 10, length.out = 1000)
#' plot(x, dgamma_mean(x, mu = 2, a = 2), type = "l")
dgamma_mean <- function(x, mu, a, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("gamma is only defined for x > 0")
  }
  if (isTRUE(a <= 0)) {
    stop("gamma is only defined for a > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("gamma is only defined for mu > 0")
  }
  return(dgamma(x = x, shape = a, rate = a / mu, log = log))
}


#' Mean parameterization of the gamma RNG.
#'
#' @param n Number of samples, scalar natural number.
#' @param mu Mean parameter, mu > 0.
#' @param a Shape parameter, a > 0.
#'
#' @return n Gamma distributed samples.
#' @export
#'
#' @examples hist(log(rgamma_mean(10000, mu = 2, a = 1)))
rgamma_mean <- function(n, mu = 1, a = 1) {
  # check the arguments
  if (isTRUE(a <= 0)) {
    stop("gamma is only defined for a > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("gamma is only defined for mu > 0")
  }
  return(rgamma(n = n, shape = a, rate = a / mu))
}
