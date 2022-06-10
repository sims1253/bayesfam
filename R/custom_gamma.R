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
dgamma_custom <- function(x, mu, a) {
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
  return(dgamma(x = x, shape = a, rate = a / mu))
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
rgamma_custom <- function(n, mu, a) {
  # check the arguments
  if (isTRUE(a <= 0)) {
    stop("gamma is only defined for a > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("gamma is only defined for mu > 0")
  }
  return(rgamma(n = n, shape = a, rate = a / mu))
}
