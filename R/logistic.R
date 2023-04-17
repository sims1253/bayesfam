#' Probability density function for the logistic distribution
#'
#' @details The logistic distribution has density
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
#' plot(x, dlogistic(x, mu = 2, sigma = 1), type = "l")
dlogistic <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(sigma <= 0))) {
    stop("logistic is only defined for sigma > 0")
  }
  # Maybe overkill?
  if (!lenEqual(list_of_vectors = list(x, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("logistic argument vectors could not be matched. May be due to wrong type,
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

#' Quantile function of the logistic distribution
#'
#' @param p quantile value, 0 < p < 1
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#'
#' @return Quantiles of the logistic distribution
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 100)
#' plot(x, qlogistic(x, mu = 2, sigma = 2), type = "l")
qlogistic <- function(p, mu, sigma) {
  # check the arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("quantile value has to be between 0 and 1")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("logistic is only defined for sigma > 0")
  }
  if (!lenEqual(list_of_vectors = list(p, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("logistic argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }

  return(mu + sigma * (log(p) - log1p(-p)))
}

#' RNG for the logistic distribution
#'
#' @param n Number of samples.
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#'
#' @return Random numbers from the logistic distribution.
#' @export
#'
#' @examples hist(rlogistic(100, mu = 2, sigma = 2))
rlogistic <- function(n, mu = 0, sigma = 1) {
  # check the arguments
  if (!isNat_len(n)) {
    stop("The number RNG-samples has to be a scalar natural")
  }
  return(qlogistic(runif(n, min = 0, max = 1), mu, sigma))
}

#' Log-Likelihood of the logistic distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of logistic given data in prep
log_lik_logistic <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dlogistic(y, mu, sigma, log = TRUE))
}


#' posterior_predict for the logistic distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ... catchall argument
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_logistic <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rlogistic(prep$ndraws, mu, sigma))
}

#' posterior_epred for the logistic distribution
#'
#' @param prep BRMS data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_logistic <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' logistic brms custom family
#'
#' @param link Link function for function
#' @param link_sigma Link function for sigma argument
#'
#' @return BRMS logistic distribution family
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = rlogistic(n = 1000, mu = a + 2, sigma = 2))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'   family = logistic(), stanvars = logistic()$stanvars,
#'   refresh = 0)
#' plot(fit)
logistic <- function(link = "identity", link_sigma = "log") {
  family <- brms::custom_family(
    "logistic_r",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(-NA, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_logistic,
    posterior_predict = posterior_predict_logistic,
    posterior_epred = posterior_epred_logistic
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real logistic_r_lpdf(real y, real mu, real sigma) {
        return logistic_lpdf(y | mu, sigma);
      }

      real logistic_r_rng(real mu, real sigma) {
        return logistic_rng(mu, sigma);
      }",
    block = "functions"
  )
  return(family)
}
