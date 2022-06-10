#' Title
#'
#' @param x
#' @param mu
#' @param sigma
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dcloglognormal <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  logpdf <- (-(log(sigma) + 0.5 * (log(2) + log(pi)))) +
    -log(-1 * ((1 - x) * log(1 - x))) +
    (-(cloglog(x) - mu)^2) / (2 * (sigma^2))
  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Title
#'
#' @param n
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
rcloglognormal <- function(n, mu, sigma) {
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    inv_cloglog(rnorm(n, mu, sigma))
  )
}

#' Title
#'
#' @param i
#' @param prep
#'
#' @return
#'
#'
#' @examples
log_lik_cloglognormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dcloglognormal(y, mu, sigma, log = TRUE))
}

#' Title
#'
#' @param i
#' @param prep
#' @param ...
#'
#' @return
#'
#'
#' @examples
posterior_predict_cloglognormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rcloglognormal(prep$ndraws, mu, sigma))
}

#' Title
#'
#' @param prep
#'
#' @return
#'
#'
#' @examples
posterior_epred_cloglognormal <- function(prep) {
  # https://doi.org/10.1080/03610926.2020.1752723 might solve this
  stop("Due to the mean not having an analytical solution for the cloglog-normal
        distribution, posterior_epred is currently not supported.")
}

#' Title
#'
#' @param link
#' @param link_sigma
#'
#' @return
#' @export
#'
#' @examples
cloglognormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "cloglognormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_cloglognormal,
    posterior_predict = posterior_predict_cloglognormal,
    posterior_epred = posterior_epred_cloglognormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real cloglognormal_lpdf(real y, real mu, real sigma) {
        return   (-(log(sigma) + 0.5 * (log(2) + log(pi())))) +
                 -log(-1 * ((1-y)*log(1-y))) +
                 (-(log(-log1m(y)) - mu)^2) / (2 * (sigma^2));
      }

      real cloglognormal_rng(real mu, real sigma) {
        return inv_cloglog(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
