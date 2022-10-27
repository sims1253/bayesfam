#' Median parameterization of the Weibull pdf.
#'
#' @param x Value space, x > 0.
#' @param mu Median parameter, mu > 0.
#' @param k Shape parameter, k > 0.
#' @param log Optional argument. If TRUE, returns log(pdf). Normally False.
#'
#' @details Define constant sigma as
#' \deqn{\sigma(\mu, k) := \mu / \Gamma(1 + 1 / k)}
#' @details The Weibull distribution density is defined as
#' \deqn{f(y) = \frac{k}{\sigma} * (\frac{x}{\sigma})^{\alpha - 1} * exp(-(\frac{x}{\sigma})^\alpha)}
#'
#' @return f(x | mu, k)
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 10, length.out = 1000)
#' plot(x, dweibull_median(x, mu = 2, k = 1), type = "l")
dweibull_median <- function(x, mu, k, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("weibull is only defined for x > 0")
  }
  if (isTRUE(k <= 0)) {
    stop("weibull is only defined for k > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("weibull is only defined for mu > 0")
  }
  return(dweibull(x = x, shape = k, scale = mu / gamma(1 + 1 / k), log))
}


#' Median parameterization of the Weibull RNG.
#'
#' @param n Number of samples, scalar natural number.
#' @param mu Median parameter, mu > 0.
#' @param k Shape parameter, k > 0.
#'
#' @return n Weibull distributed samples.
#' @export
#'
#' @examples hist(log(rweibull_mean(10000, mu = 2, k = 1)))
rweibull_median <- function(n, mu, k) {
  # check the arguments
  if (isTRUE(k <= 0)) {
    stop("weibull is only defined for k > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("weibull is only defined for mu > 0")
  }
  return(rweibull(n = n, shape = k, scale = mu / gamma(1 + 1 / k)))
}
