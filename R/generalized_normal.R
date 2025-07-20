#' Probability density function for the generalized_normal distribution
#'
#' @details The generalized_normal distribution has density
#' \deqn{f(y | \mu, \sigma, \beta) = \frac{\beta}{2 \beta \Gamma(1/\beta)}exp(-|z|^\beta)}
#' @details Where z is the linear transformation
#' \deqn{z(y, \mu, \sigma) = \frac{y - \mu}{\sigma}}
#'
#' @param x Value, unbound
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#' @param beta shape, beta > 0
#' @param log Optional argument. If true, returns the log density.
#'
#' @return Density of the pdf given x, mu and sigma
#' @export
#'
#' @examples x <- seq(from = -5, to = 10, length.out = 1000)
#' plot(x, dgeneralized_normal(x, mu = 1, sigma = 1, beta = 0.5), type = "l")
dgeneralized_normal <- function(x, mu, sigma, beta, log = FALSE) {
  # check the arguments
  if (isTRUE(any(sigma <= 0))) {
    stop("generalized_normal is only defined for sigma > 0")
  }
  if (isTRUE(any(beta <= 0))) {
    stop("generalized_normal is only defined for beta > 0")
  }

  lpdf <- log(beta) -
    (log(2) + log(sigma) + lgamma(1 / beta)) -
    (abs(x - mu) / sigma)^beta

  # return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function of the generalized_normal distribution
#'
#' @param p quantile value, 0 < p < 1
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#' @param beta shape, beta > 0
#'
#' @return Quantiles of the generalized_normal distribution
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 100)
#' plot(x, qgeneralized_normal(x, mu = 2, sigma = 2, beta = 1), type = "l")
qgeneralized_normal <- function(p, mu, sigma, beta) {
  # check the arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("quantile value has to be between 0 and 1")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("generalized_normal is only defined for sigma > 0")
  }
  if (isTRUE(any(beta <= 0))) {
    stop("generalized_normal is only defined for beta > 0")
  }

  return(
    sign(p - 0.5) *
      ((sigma^beta) *
        qgamma(2 * abs(p - 0.5), 1 / beta))^(1 / beta) +
      mu
  )
}

#' RNG for the generalized_normal distribution
#'
#' @param n Number of sample
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#' @param beta shape, beta > 0
#'
#' @return Random numbers from the generalized_normal distribution.
#' @export
#'
#' @examples hist(rgeneralized_normal(100, mu = 2, sigma = 2, beta = 2))
rgeneralized_normal <- function(n, mu = 0, sigma = 1, beta = 1) {
  return(qgeneralized_normal(
    p = runif(n, min = 0, max = 1),
    mu = mu,
    sigma = sigma,
    beta = beta
  ))
}

#' Log-Likelihood of the generalized_normal distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of generalized_normal given data in prep
log_lik_generalized_normal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dgeneralized_normal(
    x = y,
    mu = mu,
    sigma = sigma,
    beta = beta,
    log = TRUE
  ))
}


#' posterior_predict for the generalized_normal distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ... catchall argument
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_generalized_normal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  return(rgeneralized_normal(
    n = prep$ndraws,
    mu = mu,
    sigma = sigma,
    beta = beta
  ))
}

#' posterior_epred for the generalized_normal distribution
#'
#' @param prep BRMS data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_generalized_normal <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}

#' Generalized Normal BRMS family
#'
#' @param link Link function for function
#' @param link_sigma Link function for sigma argument
#' @param link_beta Link function for beta argument
#'
#' @return BRMS generalized_normal distribution family
#' @export
#'
#' @examples data <- list(y = rgeneralized_normal(n = 1000, mu = 2, sigma = 2, beta = 4))
#' fit <- brms::brm(
#'   formula = y ~ 1, data = data,
#'   family = generalized_normal(), stanvars = generalized_normal()$stanvars,
#'   init = 0.1
#' )
#' plot(fit)
generalized_normal <- function(
  link = "identity",
  link_sigma = "log",
  link_beta = "log"
) {
  family <- brms::custom_family(
    "generalized_normal",
    dpars = c("mu", "sigma", "beta"),
    links = c(link, link_sigma, link_beta),
    lb = c(-NA, 0, 0),
    ub = c(NA, NA, NA),
    type = "real",
    log_lik = log_lik_generalized_normal,
    posterior_predict = posterior_predict_generalized_normal,
    posterior_epred = posterior_epred_generalized_normal
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real generalized_normal_lpdf(real y, real mu, real sigma, real beta) {
        return log(beta) - (log(2) + log(sigma) + lgamma(1/beta))
          - (abs(y - mu) / sigma)^beta;
      }

      int sign(real x) {
        if (x > 0) {
          return 1;
        } else if (x < 0) {
          return -1;
        } else {
          return 0;
        }
      }

      real generalized_normal_rng(real mu, real sigma, real beta) {
        real p = uniform_rng(0,1);
        real q_part = ((sigma^beta) * exp(gamma_lccdf(2*abs(p - 0.5) | 1/beta, 1)))^(1/beta);
        return sign(p - 0.5) * q_part + mu;
      }",
    block = "functions"
  )
  return(family)
}
