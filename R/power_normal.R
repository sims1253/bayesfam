# Implementation after Masashi Goto and Toshiaki Inoue
# Some Properties of the Power Normal Distribution
#

#' Probability density function for the power_normal distribution
#'
#' @details The beta prime distribution has density
#' \deqn{f(y | \mu, \sigma) = \frac{e^{-z}}{\sigma(1 + e^{-z})^2}}
#' @details Where z is the linear transformation
#' \deqn{z(y, \mu, \sigma) = \frac{y - \mu}{\sigma}}
#'
#' @param x Value, unbound
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#' @param log Optional argument. If true, returns the log density.
#'
#' @return Density of the pdf given x, mu and sigma
#' @export
#'
#' @examples x <- seq(from = -5, to = 10, length.out = 1000)
#' plot(x, dpower_normal(x, mu = 2, sigma = 1), type = "l")
dpower_normal <- function(x, mu, sigma, log = FALSE) {
  stop("tbd")
  # check the arguments
  if (isTRUE(any(sigma <= 0))) {
    stop("power_normal is only defined for sigma > 0")
  }
  # Maybe overkill?
  if (!lenEqual(list_of_vectors = list(x, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("power_normal argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }
  if (!isLogic_len(log)) {
    stop("the log argument of a density has to be a scalar boolean")
  }


  z <- (x - mu) / sigma

  lpdf <- -z - log(sigma) - 2 * log1p(exp(-z))

  #return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function of the power_normal distribution
#'
#' @param p quantile value, 0 < p < 1
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#'
#' @return Quantiles of the power_normal distribution
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 100)
#' plot(x, qpower_normal(x, mu = 2, sigma = 2), type = "l")
qpower_normal <- function(p, mu, sigma) {
  # check the arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("quantile value has to be between 0 and 1")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("power_normal is only defined for sigma > 0")
  }
  if (!lenEqual(list_of_vectors = list(p, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("power_normal argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }

  q <- mu + sigma * (log(p) - log1p(-p))
  return(q)
}

#' RNG for the power_normal distribution
#'
#' @param n Number of samples.
#' @param mu Mean, mu > 0.
#' @param sigma Scale, sigma > 0
#'
#' @return Random numbers from the power_normal distribution.
#' @export
#'
#' @examples hist(rpower_normal(100, mu = 2, sigma = 2))
rpower_normal <- function(n, mu, sigma) {
  # check the arguments
  if (!isNat_len(n)) {
    stop("The number RNG-samples has to be a scalar natural")
  }
  if (!lenEqual(list_of_vectors = list(mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("power_normal argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }
  return(qpower_normal(runif(n, min = 0, max = 1), mu, sigma))
}

#' Log-Likelihood of the power_normal distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of power_normal given data in prep
log_lik_power_normal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dpower_normal(y, mu, sigma, log = TRUE))
}


#' posterior_predict for the power_normal distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_power_normal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rpower_normal(prep$ndraws, mu, sigma))
}

#' posterior_epred for the power_normal distribution
#'
#' @param prep BRMS data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_power_normal <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' power_normal brms custom family
#'
#' @param link Link function for function
#' @param link_sigma Link function for sigma argument
#'
#' @return BRMS power_normal distribution family
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = rpower_normal(n = 1000, mu = a + 2, sigma = 2))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'   family = power_normal(), stanvars = power_normal()$stanvars,
#'   refresh = 0)
#' plot(fit)
power_normal <- function(link = "identity", link_sigma = "log") {
  family <- brms::custom_family(
    "power_normal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(NA, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_power_normal,
    posterior_predict = posterior_predict_power_normal,
    posterior_epred = posterior_epred_power_normal
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real power_normal_lpdf(real y, real mu, real sigma) {
        return power_normal_lpdf(y | mu, sigma);
      }

      real power_normal_rng(real mu, real sigma) {
        return power_normalng(mu, sigma);
      }",
    block = "functions"
  )
  return(family)
}
