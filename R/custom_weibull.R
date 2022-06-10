#' Title
#'
#' @param x
#' @param mu Mean
#' @param k Shape
#'
#' @return
#' @export
#'
#' @examples
dweibull_custom <- function(x, mu, k) {
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
  return(dweibull(x = x, shape = k, scale = mu / gamma(1 + 1 / k)))
}


#' Title
#'
#' @param n
#' @param mu Mean
#' @param k Shape
#'
#' @return
#' @export
#'
#' @examples
rweibull_custom <- function(n, mu, k) {
  # check the arguments
  if (isTRUE(k <= 0)) {
    stop("weibull is only defined for k > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("weibull is only defined for mu > 0")
  }
  return(rweibull(n = n, shape = k, scale = mu / gamma(1 + 1 / k)))
}
