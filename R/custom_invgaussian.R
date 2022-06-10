#' Title
#'
#' @param x
#' @param mu
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
dinversegaussian_custom <- function(x, mu, shape, log = FALSE){
  if (isTRUE(any(x <= 0))) {
    stop("Inverse Gaussian density is only defined for x > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("Inverse Gaussian density is only defined for mu > 0")
  }
  if (isTRUE(shape <= 0)) {
    stop("Inverse Gaussian density is only defined for shape > 0")
  }
  lpdf <- 0.5 *(log(shape) - ( log(2) + log(pi) + 3 * log(x))) -
          shape * (x - mu)^2/(2 * mu^2 * x)
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Title
#'
#' @param n
#' @param mu
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
rinversegaussian_custom <- function(n, mu, shape){
  if (isTRUE(mu <= 0)) {
    stop("Inverse Gaussian density is only defined for mu > 0")
  }
  if (isTRUE(shape <= 0)) {
    stop("Inverse Gaussian density is only defined for shape > 0")
  }
  rinv_gaussian(n, mu, shape)
}


#' Title
#'
#' @param i
#' @param prep
#'
#' @return
#' @export
#'
#' @examples
log_lik_inversegaussian_custom <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  y <- prep$data$Y[i]
  return(dinversegaussian_custom(y, mu, shape, log = TRUE))
}


#' Title
#'
#' @param i
#' @param prep
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
posterior_predict_inversegaussian_custom <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  return(rgompertz(prep$ndraws, mu, shape))
}


#' Title
#'
#' @param prep
#'
#' @return
#' @export
#'
#' @examples
posterior_epred_inversegaussian_custom <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  return(mu)
}



#' Title
#'
#' @param link
#' @param link_shape
#'
#' @return
#' @export
#'
#' @examples
inversegaussian_custom <- function(link = "log", link_shape = "log") {
  family <- brms::custom_family(
    "inversegaussian_custom",
    dpars = c("mu", "shape"),
    links = c(link, link_shape),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_inversegaussian_custom,
    posterior_predict = posterior_predict_inversegaussian_custom,
    posterior_epred = posterior_epred_inversegaussian_custom
  )
  family$stanvars <- brms::stanvar(
    scode = "
        real inversegaussian_custom_lpdf(real y, real mu, real shape) {
          return (
            0.5 *(log(shape) - (log(2) + log(pi()) + 3 * log(y))) -
            shape * (y - mu)^2/(2 * mu^2 * y)
          );
        }",
    block = "functions"
  )
  return(family)
}

