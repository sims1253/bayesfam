#' Probability density function for the power_normal distribution
#' Implementation after Nelson Wayne
#' \url{https://archive.org/details/computerprogramp4760nels/page/10/mode/2up}
#' @details The beta prime distribution has density
#' \deqn{f(y | \mu, \sigma) = \frac{e^{-z}}{\sigma(1 + e^{-z})^2}}
#' @details Where z is the linear transformation
#' \deqn{z(y, \mu, \sigma) = \frac{y - \mu}{\sigma}}
#'
#' @param x Value, unbound
#' @param mu Median, unbound
#' @param sigma Scale, sigma > 0
#' @param beta Shape, beta > 0
#' @param log Optional argument. If true, returns the log density.
#'
#' @return Density of the pdf given x, mu and sigma
#' @export
#'
#' @examples x <- seq(from = -4, to = 4, length.out = 1000)
#' plot(x, dpower_normal(x, mu = 0, sigma = 1, beta = 4), type = "l")
dpower_normal <- function(x, mu, sigma, beta, log = FALSE) {
  # check the arguments
  if (isTRUE(any(sigma <= 0))) {
    stop("power_normal is only defined for sigma > 0")
  }
  if (isTRUE(any(beta <= 0))) {
    stop("generalized_normal is only defined for beta > 0")
  }
  # Maybe overkill?
  if (!lenEqual(list_of_vectors = list(x, mu, sigma, beta), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("power_normal argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }
  if (!isLogic_len(log)) {
    stop("the log argument of a density has to be a scalar boolean")
  }

  #loc <- mu - power_zf(0.5, beta) * sigma
  #z <- (x - loc) / sigma
  z <- (x - mu)/sigma + power_zf(0.5, beta)

  # optimized version of log cdf?
  lpdf <- log(beta) - log(sigma) + dnorm(z, log = TRUE) + (beta-1) * log(pnorm(-z, lower.tail = FALSE))

  #return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Power transformed qnorm (zf-parameter used by Nelson Wayne)
#'
#' @param p percentile to be calculated, 0 < p < 1
#' @param beta shape argument, beta > 0
#'
#' @return zf argument used in power-normal computation
#'
#' @examples
power_zf <- function(p, beta) {
  # only gets called by internal functions, so no further checks necassary!
  arg <- (1 - p) ^ (1 / beta)
  return(qnorm(-arg + 1))
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
#' plot(x, qpower_normal(x, mu = 0, sigma = 2, beta = 10), type = "l")
qpower_normal <- function(p, mu, sigma, beta) {
  # check the arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("quantile value has to be between 0 and 1")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("power_normal is only defined for sigma > 0")
  }
  if (!lenEqual(list_of_vectors = list(p, mu, sigma, beta), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("power_normal argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }

  # loc <- mu - power_zf(0.5, beta) * sigma
  # q <- loc + power_zf(p, beta) * sigma

  # loc <- mu - power_zf(0.5, beta) * sigma
  q <- mu + sigma * (power_zf(p, beta) - power_zf(0.5, beta))
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
#' @examples hist(rpower_normal(100, mu = 0, sigma = 2, beta = 10))
rpower_normal <- function(n, mu, sigma, beta) {
  # check the arguments
  if (!isNat_len(n)) {
    stop("The number RNG-samples has to be a scalar natural")
  }
  return(qpower_normal(runif(n, min = 0, max = 1), mu, sigma, beta))
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
  beta <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dpower_normal(y, mu, sigma, beta, log = TRUE))
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
  beta <- brms::get_dpar(prep, "beta", i = i)
  return(rpower_normal(prep$ndraws, mu, sigma, beta))
}

#' posterior_epred for the power_normal distribution
#'
#' @param prep BRMS data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_power_normal <- function(prep) {
  #return(brms::get_dpar(prep, "mu"))
  stop("No closed form for mean")
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
#' data <- list(a = a, y = rpower_normal(n = 1000, mu = a + 1, sigma = 2, beta = 10))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'   family = power_normal(), stanvars = power_normal()$stanvars,
#'   refresh = 0)
#' plot(fit)
power_normal <- function(link = "identity", link_sigma = "log", link_beta = "log") {
  family <- brms::custom_family(
    "power_normal",
    dpars = c("mu", "sigma", "beta"),
    links = c(link, link_sigma, link_beta),
    lb = c(NA, 0, 0),
    ub = c(NA, NA, NA),
    type = "real",
    log_lik = log_lik_power_normal,
    posterior_predict = posterior_predict_power_normal,
    posterior_epred = posterior_epred_power_normal
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real power_zf(real p, real beta) {
        real arg = (1 - p) ^ (1 / beta);
        return exp(normal_lccdf(-arg + 1 | 0, 1));
      }

      real power_normal_lpdf(real y, real mu, real sigma, real beta) {
        //real loc = mu - power_zf(0.5, beta) * sigma;
        //real z = (y - loc) / sigma;

        real z = (y - mu) / sigma + power_zf(0.5, beta);
        return log(beta) - log(sigma) + normal_lpdf(z | 0, 1) + (beta-1) * normal_lcdf(-z | 0, 1);
      }

      real power_normal_rng(real mu, real sigma, real beta) {
        //real loc = mu - power_zf(0.5, beta) * sigma;
        //return loc + power_zf(p, beta) * sigma;

        real p = uniform_rng(0,1);
        return mu + sigma * (power_zf(p, beta) - power_zf(0.5, beta));
      }",
    block = "functions"
  )
  return(family)
}
