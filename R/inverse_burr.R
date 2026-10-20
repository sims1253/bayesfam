#' Probability density function for the Inverse Burr distribution, with Median parametrization.
#'
#' @source Parameterization follows the actuar package (also known as the
#' Dagum distribution), with \eqn{\tau} = shape1, \eqn{\gamma} = shape2 and
#' scale \eqn{\sigma = \mu (2^{1/\tau} - 1)^{1/\gamma}}.
#'
#' @param x Model space, defined for x > 0
#' @param mu Median parameter of pdf, mu > 0
#' @param tau Shape1 parameter of pdf, tau > 0
#' @param gamma Shape2 parameter of pdf, gamma > 0
#' @param log Optional log argument, if true, return log(pdf)
#'
#' @details \deqn{f(y | \mu, \tau, \gamma) = \frac{\tau \gamma (y/\sigma)^{\gamma \tau}}{y [1 + (y/\sigma)^{\gamma}]^{\tau + 1}}}
#'
#' @return PDF of Inverse Burr Distribution
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 10, length.out = 100)
#' plot(x, dinverse_burr(x, mu = 2, tau = 2, gamma = 3), type = "l")
dinverse_burr <- function(x, mu, tau, gamma, log = FALSE) {
  # check arguments
  if (isTRUE(any(x <= 0))) {
    stop("The Inverse Burr PDF is defined only on the positive scale")
  }
  if (isTRUE(any(mu <= 0 | tau <= 0 | gamma <= 0))) {
    stop("The Inverse Burr PDF is defined only for mu, tau, gamma > 0")
  }
  scale <- inverse_burr_scale(mu, tau, gamma)
  lpdf <- log(tau) +
    log(gamma) +
    tau * gamma * (log(x) - log(scale)) -
    log(x) -
    (tau + 1) * log1p((x / scale)^gamma)

  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Cumulative distribution function for the Inverse Burr distribution, with Median parametrization.
#'
#' @param q Value to evaluate the CDF at, q > 0
#' @param mu Median parameter, mu > 0
#' @param tau Shape1 parameter, tau > 0
#' @param gamma Shape2 parameter, gamma > 0
#'
#' @return CDF of Inverse Burr Distribution
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 10, length.out = 100)
#' plot(x, pinverse_burr(x, mu = 2, tau = 2, gamma = 3), type = "l")
pinverse_burr <- function(q, mu, tau, gamma) {
  # check arguments
  if (isTRUE(any(q <= 0))) {
    stop("The Inverse Burr CDF is defined only on the positive scale")
  }
  if (isTRUE(any(mu <= 0 | tau <= 0 | gamma <= 0))) {
    stop("The Inverse Burr CDF is defined only for mu, tau, gamma > 0")
  }
  w <- (q / inverse_burr_scale(mu, tau, gamma))^gamma
  return((w / (1 + w))^tau)
}

#' Quantile function for the Inverse Burr distribution, with Median parametrization.
#'
#' @param p Quantile to be calculated, p in (0, 1)
#' @param mu Median parameter, mu > 0
#' @param tau Shape1 parameter, tau > 0
#' @param gamma Shape2 parameter, gamma > 0
#'
#' @return Inverse of CDF, calculates a value, given a probability p
#' @export
#'
#' @examples p <- seq(from = 0.01, to = 0.99, length.out = 100)
#' plot(p, qinverse_burr(p, mu = 2, tau = 2, gamma = 3), type = "l")
qinverse_burr <- function(p, mu, tau, gamma) {
  # check arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("The Inverse Burr quantile function is defined only for p in (0, 1)")
  }
  if (isTRUE(any(mu <= 0 | tau <= 0 | gamma <= 0))) {
    stop(
      "The Inverse Burr quantile function is defined only for mu, tau, gamma > 0"
    )
  }
  z <- p^(1 / tau)
  return(inverse_burr_scale(mu, tau, gamma) * (z / (1 - z))^(1 / gamma))
}

#' RNG function for the Inverse Burr distribution, with Median parametrization.
#'
#' @param n Number of draws
#' @param mu Median parameter, mu > 0
#' @param tau Shape1 parameter, tau > 0
#' @param gamma Shape2 parameter, gamma > 0
#'
#' @return An Inverse Burr distributed RNG vector of size n
#' @export
#'
#' @examples hist(log(rinverse_burr(1000, mu = 2, tau = 2, gamma = 3)))
rinverse_burr <- function(n, mu = 1, tau = 1, gamma = 1) {
  return(qinverse_burr(runif(n), mu = mu, tau = tau, gamma = gamma))
}

# actuar scale parameter for the median parametrization
#' @noRd
inverse_burr_scale <- function(mu, tau, gamma) {
  return(mu * (2^(1 / tau) - 1)^(1 / gamma))
}

#' Log-Likelihood vignette for the Inverse Burr distribution, with Median parametrization.
#'
#' @param i brms indices
#' @param prep brms data
#'
#' @return Log-Likelihood of Inverse Burr given data in prep
log_lik_inverse_burr <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  gamma <- brms::get_dpar(prep, "gamma", i = i)
  y <- prep$data$Y[i]
  return(dinverse_burr(y, mu, tau, gamma, log = TRUE))
}

#' Posterior-Prediction vignette for the Inverse Burr distribution, with Median parametrization.
#'
#' @param i brms indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Posterior prediction of Inverse Burr, given data in prep
posterior_predict_inverse_burr <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  gamma <- brms::get_dpar(prep, "gamma", i = i)
  return(rinverse_burr(prep$ndraws, mu, tau, gamma))
}

#' Expectation-Predict vignette for the Inverse Burr distribution, with Median parametrization.
#'
#' @param prep brms data
#'
#' @return Mean of the Inverse Burr distribution given data in prep; Inf for
#'   draws with gamma <= 1, where the mean does not exist
posterior_epred_inverse_burr <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  tau <- brms::get_dpar(prep, "tau")
  gamma <- brms::get_dpar(prep, "gamma")
  epred <- inverse_burr_scale(mu, tau, gamma) *
    exp(lgamma(tau + 1 / gamma) + lgamma(1 - 1 / gamma) - lgamma(tau))
  epred[gamma <= 1] <- Inf
  return(epred)
}

#' Inverse Burr Stan-implementation in Median parametrization.
#'
#' @param link Link function for function
#' @param link_tau Link function for tau argument
#' @param link_gamma Link function for gamma argument
#'
#' @return brms Inverse Burr distribution family
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(a = a, y = rinverse_burr(1000, exp(0.5 * a + 1), 2, 3))
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = inverse_burr(), stanvars = inverse_burr()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
inverse_burr <- function(link = "log", link_tau = "log", link_gamma = "log") {
  family <- brms::custom_family(
    "inverse_burr",
    dpars = c("mu", "tau", "gamma"),
    links = c(link, link_tau, link_gamma),
    lb = c(0, 0, 0),
    ub = c(NA, NA, NA),
    type = "real",
    log_lik = log_lik_inverse_burr,
    posterior_predict = posterior_predict_inverse_burr,
    posterior_epred = posterior_epred_inverse_burr
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real inverse_burr_lpdf(real y, real mu, real tau, real gamma) {
        real scale = mu * pow(pow(2, inv(tau)) - 1, inv(gamma));
        return log(tau) + log(gamma) +
               tau * gamma * (log(y) - log(scale)) -
               log(y) -
               (tau + 1) * log1p(pow(y / scale, gamma));
      }

      real inverse_burr_rng(real mu, real tau, real gamma) {
        real scale = mu * pow(pow(2, inv(tau)) - 1, inv(gamma));
        real z = pow(uniform_rng(0, 1), inv(tau));
        return scale * pow(z / (1 - z), inv(gamma));
      }",
    block = "functions"
  )
  return(family)
}
