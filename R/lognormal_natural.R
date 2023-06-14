
#' Lognormal Natural Density function with value space parameterisation (mu same as x)
#'
#' @param x Value space, x > 0
#' @param mu Mean parameter, mu > 0
#' @param sigma Shape parameter, sigma unbound
#' @param log boolean value to return log-pdf, default FALSE
#'
#' @return Density of the Lognormal Natural given x, mu and sigma.
#' @export
#'
#' @examples x <- seq(from=0.1, to=5, length.out=100)
#' plot(x, dlognormal_natural(x, 1, 2))
dlognormal_natural <- function(x, mu, sigma, log = FALSE) {
  if(isTRUE(any(mu <= 0))) {
    stop("Mu has to be > 0")
  }
  common_term <- log1p(sigma^2/mu^2)
  return(dlognormal(x, log(mu)-common_term/2, sqrt(common_term), log))
}

#' Lognormal Natural RNG function
#'
#' @param n N number of samples to draw, n is a scalar natural number
#' @param mu mean argument, mu > 0
#' @param sigma shape argument, sigma > 0
#'
#' @return N samples distributed as Lognormal Natural Likelihood
#' @export
#'
#' @examples hist(rlognormal_natural(100, 1, 2))
rlognormal_natural <- function(n, mu = 0, sigma = 1) {
  if(isTRUE(any(mu <= 0))) {
    stop("Mu has to be > 0")
  }
  common_term <- log1p(sigma^2/mu^2)
  return(rlognormal(n, log(mu)-common_term/2, sqrt(common_term)))
}

#' Logarithmic Density BRMS vignette of Lognormal Natural Likelihood
#'
#' @param i Indices of BRMS data
#' @param prep BRMS data
#'
#' @return Log-Likelihood of BRMS data
#' @export
log_lik_lognormal_natural <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)

  y <- prep$data$Y[i]
  return(dlognormal_natural(y, mu, sigma, log=TRUE))
}


#' Posterior Prediction BRMS Vignette of Lognormal Natural Likelihood
#'
#' @param i Indices of BRMS data
#' @param prep BRMS data
#' @param ... Catchall argument
#'
#' @return  Draws from the Posterior Predictive Distribution
#' @export
posterior_predict_lognormal_natural <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  return(rlognormal_natural(prep$ndraws, mu, sigma))
}

#' Posterior mean BRMS vingette of Lognormal Natural Likelihood
#'
#' @param prep BRMS data
#'
#' @return Mean of the Lognormal natural posterior
#' @export
posterior_epred_lognormal_natural <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  return(mu)
}

#' Lognormal Natural BRMS family
#'
#' @param link link for the mean argument, default = log
#' @param link_sigma link for the shape argument, default = log
#'
#' @return BRMS model object
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = rlognormal_natural(n = 1000, mu = exp(0.5 * a + 1), sigma = exp(2)))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'  family = lognormal_natural(), stanvars = lognormal_natural()$stanvars,
#'  refresh = 0)
#' plot(fit)
lognormal_natural <- function(link = "log", link_sigma = "log") {
  family <- brms::custom_family(
    "lognormal_natural",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
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
