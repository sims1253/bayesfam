#' Custom rstudent with default df value
#' mu and df arguments switched to comply to the Bayesfam allowing compatibility
#' for Bayesim
#'
#' @param n Number of sampels to draw, has to be a scalar natural
#' @param mu Mean argument, mu unbound
#' @param df Degrees of freedom variable
#' @param sigma Shape parameter, sigma > 0
#'
#' @return Vector of length n in student distribution
#' @export
#'
#' @examples hist(rstudent_mean(100, 1, 2, 1))
rstudent_mean <- function(n, mu = 0, df = 1, sigma = 1) {
  if (isTRUE(any(df <= 0))) {
    stop("df has to be bigger than 0")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("sigma has to be bigger than 0")
  }

  return(brms::rstudent_t(n, df, mu, sigma))
}
