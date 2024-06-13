#' Probability density function for the discrete Weibull distribution
#'
#' @source Scholz, F., & Works, B. P. (1996). Maximum likelihood estimation for
#'         type I censored Weibull data including covariates. ISSTECH-96-022,
#'         Boeing Information and Support Services.
#'
#' @details The location-scale parameterization has the form
#' \deqn{f(x) = 1 - exp[-exp[\frac{x - \mu}{\sigma}}
#' @details Where \deqn{\mu = log(\alpha)} and \deqn{\sigma = \frac{1}{\beta}}
#'          of for the scale \deqn{\alpha} and shape \degn{\beta} of the standard
#'          parameterization:
#'          \deqn{f(x) = 1 - exp[-(\frac{x + 1}{\alpha})^\beta}
#'
#' @param x Value, x > 0.
#' @param mu Location
#' @param sigma Scale, sigma > 0.
#' @param log Optional argument. If true, returns the log density.
#'
#' @return Density of the pdf given x, mu and sigma.
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 20, length.out = 1000)
#' plot(x, dbetaprime(x, mu = 4, phi = 2), type = "l")
ddiscrete_weibull <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("The discrete Weibull distribution is only defined for x > 0")
  }
  if (!isTrue(all.equal(x, as.integer(x)))) {
    stop("The discrete Weibull distribution expects discrete inputs")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("The discrete Weibull distribution is only defined for positive sigma values")
  }

  pdf <- 1 - exp(-exp((x - mu)/sigma))

  # return either the log or the pdf itself, given the log-value
  if (log) {
    return(log1p(-pdf))
  } else {
    return(pdf)
  }
}

#' Quantile function of the discrete Weibull distribution
#'
#' @param p Probabilities, for which to calculate the quantiles
#' @param mu Location
#' @param sigma Scale, sigma > 0.
#'
#' @return Quantiles of the discrete Weibull distribution
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' plot(x, qbetaprime(x, mu = 1, phi = 2), type = "l")
qdiscrete_weibull <- function(p, mu, sigma) {
  # check the arguments
  if (isTRUE(any(p < 0 | p > 1))) {
    stop("p has to be in an interval of [0, 1]")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("The discrete Weibull distribution is only defined for positive sigma values")
  }

  q <- ceiling(-1 + exp(mu) * (-log(1 - p))^(sigma))

  return(q)
}

#' RNG for the discrete Weibull distribution
#'
#' @param n Number of samples.
#' @param mu Location
#' @param sigma Scale, sigma > 0.
#'
#' @return Random numbers from the discrete Weibull distribution.
#' @export
#'
#' @examples hist(rbetaprime(100, mu = 1, phi = 2))
rdiscrete_Weibull <- function(n, mu = 1, sigma = 1) {
  # check the arguments
  if (isTRUE(any(sigma <= 0))) {
    stop("The discrete Weibull distribution is only defined for positive sigma values")
  }
  return(qdiscrete_weibull(runif(n, min = 0, max = 1), mu, sigma))
}

#' Log-Likelihood of the discrete Weibull distribution
#'
#' @param i brms indices
#' @param prep brms data
#'
#' @return Log-Likelihood of discrete Weibull given data in prep
log_lik_discrete_weibull <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(ddiscrete_weibull(y, mu, sigma, log = TRUE))
}


#' posterior_predict for the discrete Weibull distribution
#'
#' @param i brms indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_discrete_weibull <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rdiscrete_Weibull(prep$ndraws, mu, sigma))
}

#' posterior_epred for the discrete Weibull distribution
#'
#' @param prep brms data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_discrete_weibull <- function(prep) {
  stop("Posterior_epred is not defined for the discrete Weibull distribution!")
}


#' Discrete Weibull brms custom family
#'
#' @param link Link function for the location parameter
#' @param link_sigma Link function for the scale parameter
#'
#' @return brms discrete Weibull distribution family
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = rbetaprime(n = 1000, mu = exp(0.5 * a + 1), phi = 2))
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = betaprime(), stanvars = betaprime()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
discrete_weibull <- function(link = "identity", link_sigma = "log") {
  family <- brms::custom_family(
    "discrete_weibull",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(NA, 0),
    ub = c(NA, NA),
    type = "int",
    log_lik = log_lik_discrete_weibull,
    posterior_predict = posterior_predict_discrete_weibull,
    posterior_epred = posterior_epred_discrete_weibull
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real discrete_weibull_lpmf(int y, real mu, real sigma) {
      return log1m(exp(-exp((y - mu)/sigma)));
      }",
    block = "functions"
  )
  return(family)
}
