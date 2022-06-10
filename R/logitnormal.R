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
dlogitnormal <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  logpdf <- (-(log(sigma) + 0.5 * (log(2) + log(pi)))) +
    (-(log(x) + log1p(-x))) +
    (-(logit(x) - mu)^2) / (2 * (sigma^2))
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
rlogitnormal <- function(n, mu, sigma) {
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    inv_logit(rnorm(n, mu, sigma))
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
log_lik_logitnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dlogitnormal(y, mu, sigma, log = TRUE))
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
posterior_predict_logitnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rlogitnormal(prep$ndraws, mu, sigma))
}

#' Title
#'
#' @param prep
#'
#' @return
#'
#'
#' @examples
posterior_epred_logitnormal <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(exp(mu + sigma^2 / 2))
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
logitnormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "logitnormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_logitnormal,
    posterior_predict = posterior_predict_logitnormal,
    posterior_epred = posterior_epred_logitnormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real logitnormal_lpdf(real y, real mu, real sigma) {
        return (-(log(sigma) + 0.5 * (log(2) + log(pi())))) +
               (-(log(y) + log1m(y))) +
               (-(logistic_lccdf(y| 0, 1) - mu)^2) / (2 * (sigma^2));
      }

      real logitnormal_rng(real mu, real sigma) {
        return inv_logit(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
