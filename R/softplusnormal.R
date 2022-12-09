#' Softplus density distribution in median parametrization.
#'
#' @param x Value space of the distribution, x > 0
#' @param mu Median parameter, mu is already log-transformed, mu unbound
#' @param sigma Sigma shape parameter, sigma >= 0
#' @param log Bool argument, if true, returns the logarithmic density
#'
#' @return Normal distribution density with logit link function
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 10, length.out = 1000)
#' plot(x, dsoftplusnormal(x, mu = 1, sigma = 2), type = "l")
dsoftplusnormal <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("softplusnormal is only defined for x > 0")
  }
  if (isTRUE(sigma <= 0)) {
    stop("softplusnormal is only defined for sigma > 0")
  }
  logpdf <-
    -(log(sigma) + 0.5 * log(2 * pi)) +
    x - log(exp(x) - 1) +
    -0.5 * ((log(exp(x) - 1) - mu) / sigma)^2
  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Softplus RNG-function in median parametrization.
#'
#' @param n Number of draws
#' @param mu Median parameter, mu unbound, mu already log transformed
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @returns n Softplus distributed samples
#'
#' @export
#'
#' @examples hist(rsoftplusnormal(100, 1, 2))
rsoftplusnormal <- function(n, mu, sigma) {
  # check the arguments
  if (isTRUE(sigma <= 0)) {
    stop("softplusnormal is only defined for sigma > 0")
  }
  return(
    log(exp(rnorm(n, mu, sigma)) + 1)
  )
}

#' Log-Likelihood vignette for the Softplus distribution, in Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#'
#' @return log_likelihood of the Softplus distribution, given some BRMS data.
log_lik_softplusnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dsoftplusnormal(y, mu, sigma, log = TRUE))
}

#' Posterior-predict vignette for the Softplus distribution, with Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#' @param ...
#'
#' @return The posterior prediction of the Softplus distribution, given some BRMS data.
posterior_predict_softplusnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rsoftplusnormal(prep$ndraws, mu, sigma))
}

#' Posterior expected value prediction vignette for Softplus distribution.
#'
#' @param prep BRMS data
#'
#' @return Mean of Posterior
posterior_epred_softplusnormal <- function(prep) {
  #   mu <- brms::get_dpar(prep, "mu")
  #   sigma <- brms::get_dpar(prep, "sigma")
  #   return(-0.5 * erf((0.707107 * (mu - log(1 + exp(x)))) / sigma))
  # }
  stop("No implementation of posterior_epred for the softplus-normal
       available!")
}

#' Custom BRMS family Softplus in median parametrization.
#'
#' @param link Link function argument (as string) for Median argument. Left as identity!
#' @param link_sigma Link function argument (as string) for Shape argument
#'
#' @return Softplus BRMS model-object
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(a = a, y = rsoftplusnormal(1000, 0.5 * a + 1, 2))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'  family = softplusnormal(), stanvars = softplusnormal()$stanvars,
#'  refresh = 0)
#' plot(fit)
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
      return -(log(sigma) + 0.5 * log(2 * pi())) +
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
