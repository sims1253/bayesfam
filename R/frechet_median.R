#' Median parameterization of the Fréchet pdf.
#'
#' @param x x value space, x > 0
#' @param mu Median
#' @param nu Shape, nu > 0
#'
#' @details Define scale parameter sigma as
#' \deqn{\sigma(\mu, \nu) := \mu \cdot \log(2)^{1 / \nu}}
#' @details The Frechet distribution has density
#' \deqn{f(y) = (\nu /\sigma) * (y / \sigma)^{-(\nu + 1)} * exp(-(y / \sigma)^{-\nu}) }
#'
#' @return dfrechet(x | mu, nu)
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 20, length.out = 1000)
#' plot(x, dfrechet_median(x, mu = 6, nu = 4), type = "l")
dfrechet_median <- function(x, mu, nu) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("frechet is only defined for x > 0")
  }
  if (isTRUE(nu <= 0)) {
    stop("frechet is only defined for nu > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::dfrechet(
    x = x,
    loc = 0,
    scale = mu * log(2)^(1 / nu),
    shape = nu
  ))
}


#' Median parameterization of the Fréchet RNG
#'
#' @param n Number samples to draw
#' @param mu Median
#' @param nu Shape, nu > 0
#'
#' @return n samples in Frechet-Distribution
#' @export
#'
#' @examples hist(rfrechet_median(100, mu = 1, nu = 2))
rfrechet_median <- function(n, mu = 1, nu = 2) {
  # check the arguments
  if (isTRUE(nu <= 0)) {
    stop("frechet is only defined for nu > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::rfrechet(n = n, scale = mu * log(2)^(1 / nu), shape = nu))
}
