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
dlognormal <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("lognormal is only defined for x > 0")
  }
  if (isTRUE(sigma <= 0)) {
    stop("lognormal is only defined for sigma > 0")
  }
  logpdf <-
    -(log(x) + log(sigma) + 0.5 * (log(2) + log(pi))) +
    (-(log(x) - mu)^2 / (2 * sigma^2))
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
rlognormal <- function(n, mu, sigma) {
  # check the arguments
  if (isTRUE(sigma <= 0)) {
    stop("lognormal is only defined for sigma > 0")
  }
  return(
    exp(rnorm(n, mu, sigma))
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
log_lik_lognormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dlognormal(y, mu, sigma, log = TRUE))
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
posterior_predict_lognormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rlognormal(prep$ndraws, mu, sigma))
}

#' Title
#'
#' @param prep
#'
#' @return
#'
#'
#' @examples
posterior_epred_lognormal <- function(prep) {
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
lognormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "lognormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_lognormal,
    posterior_predict = posterior_predict_lognormal,
    posterior_epred = posterior_epred_lognormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real lognormal_lpdf(real y, real mu, real sigma) {
        return -(log(y) + log(sigma) + 0.5 * (log(2) + log(pi()))) +
                (-(log(y) - mu)^2 / (2 * sigma^2));
      }

      real lognormal_rng(real mu, real sigma) {
        return exp(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
