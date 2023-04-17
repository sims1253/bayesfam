#' Custom rexgauss with default  values
#'
#' @param n Number of sampels to draw, has to be a scalar natural
#' @param mu Mean argument, mu unbound
#' @param sigma Shape parameter, sigma > 0
#' @param beta Shape parameter, beta > 0
#'
#' @return Vector of length n in exgaussian distribution
#' @export
#'
#' @examples hist(rstudent_custom(100, 1, 2, 1))
rexgauss_custom <- function(n, mu = 0, sigma = 1, beta = 1) {
  if(!isNat_len(n)) {
    stop("n has to be a scalar natural number")
  }
  if(isTRUE(any(sigma <= 0))) {
    stop("sigma has to be bigger than 0")
  }
  if(isTRUE(any(beta <= 0))) {
    stop("beta has to be bigger than 0")
  }
  if(!lenEqual(list(mu, sigma, beta), scalars_allowed = TRUE,
               type_check = is.numeric, na_allowed = FALSE)) {
    stop("Type error, or unequal length vectors with n > 1")
  }
  return(brms::rexgaussian(n,  mu, sigma, beta))
}
