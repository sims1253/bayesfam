#' Title
#'
#' @param x
#' @param mu Mean
#' @param nu Shape
#'
#' @return
#' @export
#'
#' @examples
dfrechet_custom <- function(x, mu, nu) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("frechet is only defined for x > 0")
  }
  if (isTRUE(nu <= 1)) {
    stop("frechet is only defined for nu > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::dfrechet(x = x, loc = 0, scale = mu / gamma(1 - 1 / nu), shape = nu))
}


#' Title
#'
#' @param n
#' @param mu Mean
#' @param nu Shape
#'
#' @return
#' @export
#'
#' @examples
rfrechet_custom <- function(n, mu, nu) {
  # check the arguments
  if (isTRUE(nu <= 1)) {
    stop("frechet is only defined for nu > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("frechet is only defined for mu > 0")
  }
  return(brms::rfrechet(n = n, scale = mu / gamma(1 - 1 / nu), shape = nu))
}
