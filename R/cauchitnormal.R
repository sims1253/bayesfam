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
dcauchitnormal <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  logpdf <- (-(log(sigma) + 0.5 * (log(2) + log(pi)))) +
    log(pi) + 2 * (-log(cos(pi * (x - 0.5)))) +
    (-(cauchit(x) - mu)^2) / (2 * (sigma^2))
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
rcauchitnormal <- function(n, mu, sigma) {
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    inv_cauchit(rnorm(n, mu, sigma))
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
log_lik_cauchitnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dcauchitnormal(y, mu, sigma, log = TRUE))
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
posterior_predict_cauchitnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rcauchitnormal(prep$ndraws, mu, sigma))
}

#' Title
#'
#' @param prep
#'
#' @return
#'
#'
#' @examples
posterior_epred_cauchitnormal <- function(prep) {
  # https://doi.org/10.1080/03610926.2020.1752723 might solve this
  stop("Due to the mean not having an analytical solution for the cauchit-normal
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
cauchitnormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "cauchitnormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_cauchitnormal,
    posterior_predict = posterior_predict_cauchitnormal,
    posterior_epred = posterior_epred_cauchitnormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real cauchitnormal_lpdf(real y, real mu, real sigma) {
        return (-(log(sigma) + 0.5 * (log(2) + log(pi())))) +
               log(pi()) + 2 * (-log(cos(pi() * (y - 0.5)))) +
               (-(cauchy_lccdf(y| 0, 1) - mu)^2) / (2 * (sigma^2));
      }

      real cauchitnormal_rng(real mu, real sigma) {
        return cauchy_cdf(normal_rng(mu, sigma), 0, 1);
      }",
    block = "functions"
  )
  return(family)
}
