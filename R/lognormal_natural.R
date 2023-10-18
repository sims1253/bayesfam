#' Lognormal Natural density function with value space parameterisation (mu same as x)
#'
#' @param x vector of values, x > 0
#' @param mu Mean, mu > 0
#' @param sigma Shape, sigma > 0
#' @param log logical; if TRUE, log(pdf) is returned
#'
#' @return Density of the Lognormal Natural given x, mu and sigma.
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 5, length.out = 100)
#' plot(x, dlognormal_natural(x, 1, 2))
dlognormal_natural <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(mu <= 0))) {
    stop("Mu has to be > 0")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("sigma has to be > 0")
  }
  common_term <- log1p(sigma^2 / mu^2)
  return(dlognormal(x, log(mu) - common_term / 2, sqrt(common_term), log))
}

#' Lognormal Natural RNG function
#'
#' @param n number of observations
#' @param mu mean, mu > 0
#' @param sigma sigma, sigma > 0
#'
#' @return n samples drawn from the Lognormal natural distribution
#' @export
#'
#' @examples hist(rlognormal_natural(100, 1, 2))
rlognormal_natural <- function(n, mu = 1, sigma = 1) {
  if (isTRUE(any(mu <= 0))) {
    stop("Mu has to be > 0")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("sigma has to be > 0")
  }
  common_term <- log1p(sigma^2 / mu^2)
  return(rlognormal(n, log(mu) - common_term / 2, sqrt(common_term)))
}

#' Log-Likelihood of the Lognormal Natural family
#'
#' @param i Indices
#' @param prep brms data
#'
#' @return Log-Likelihood of brms data
log_lik_lognormal_natural <- function(i, prep) {
  return(
    dlognormal_natural(
      x = prep$data$Y[i],
      mu = brms::get_dpar(prep, "mu", i = i),
      sigma = brms::get_dpar(prep, "sigma", i = i),
      log = TRUE
    )
  )
}


#' Posterior Prediction of the Lognormal Natural family
#'
#' @param i Indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Draws from the posterior predictive distribution
posterior_predict_lognormal_natural <- function(i, prep, ...) {
  return(
    rlognormal_natural(
      n = prep$ndraws,
      mu = brms::get_dpar(prep, "mu", i = i),
      sigma = brms::get_dpar(prep, "sigma", i = i)
    )
  )
}

#' Expectations of posterior predictions of the Lognormal Natural family
#'
#' @param prep brms data
#'
#' @return Mean of the posterior predictive distribution
posterior_epred_lognormal_natural <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}

#' Lognormal Natural brms family
#'
#' @param link link for mu, default = log
#' @param link_sigma link for sigma, default = log
#'
#' @return lognormal natural brms family
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = rlognormal_natural(n = 1000, mu = exp(0.5 * a + 1), sigma = exp(2)))
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = lognormal_natural(), stanvars = lognormal_natural()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
lognormal_natural <- function(link = "log", link_sigma = "log") {
  family <- brms::custom_family(
    "lognormal_natural",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_lognormal_natural,
    posterior_predict = posterior_predict_lognormal_natural,
    posterior_epred = posterior_epred_lognormal_natural
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real lognormal_natural_lpdf(real y, real mu, real sigma) {
        real common_term = log1p(sigma^2/mu^2);
        return lognormal_lpdf(y | log(mu)-common_term/2,
                                  sqrt(common_term));
      }
      real lognormal_natural_rng(real mu, real sigma) {
        real common_term = log1p(sigma^2/mu^2);
        return lognormal_rng(log(mu)-common_term/2,
                                sqrt(common_term));
      }
    ",
    block = "functions"
  )
  return(family)
}
