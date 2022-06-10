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
dsoftplusnormal <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("softplusnormal is only defined for x > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("softplusnormal is only defined for mu > 0")
  }
  if (isTRUE(sigma <= 0)) {
    stop("softplusnormal is only defined for sigma > 0")
  }
  logpdf <-
    -(log(sigma) + 0.5 * (log(2) - log(pi))) +
    x - log(exp(x) - 1) +
    -0.5 * ((log(exp(x) - 1) - mu) / sigma)^2
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
rsoftplusnormal <- function(n, mu, sigma) {
  # check the arguments
  if (isTRUE(mu <= 0)) {
    stop("softplusnormal is only defined for mu > 0")
  }
  if (isTRUE(sigma <= 0)) {
    stop("softplusnormal is only defined for sigma > 0")
  }
  return(
    log(exp(rnorm(n, mu, sigma)) + 1)
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
log_lik_softplusnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dsoftplusnormal(y, mu, sigma, log = TRUE))
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
posterior_predict_softplusnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rsoftplusnormal(prep$ndraws, mu, sigma))
}

#' Title
#'
#' @param prep
#'
#' @return
#'
#'
#' @examples
posterior_epred_softplusnormal <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(-0.5 * erf((0.707107 * (mu - log(1 + exp(x)))) / sigma))
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
softplusnormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "softplusnormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_softplusnormal,
    posterior_predict = posterior_predict_softplusnormal,
    posterior_epred = posterior_epred_softplusnormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real softplusnormal_lpdf(real y, real mu, real sigma) {
      return -(log(sigma) + 0.5 * (log(2) - log(pi()))) +
              y - log(exp(y) - 1) +
              -0.5 * ((log(exp(y) - 1) - mu)/sigma)^2;
      }

      real softplusnormal_rng(real mu, real sigma) {
        return log(exp(normal_rng(mu, sigma)) - 1);
      }",
    block = "functions"
  )
  return(family)
}
