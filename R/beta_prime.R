#' Probability density function for the beta prime distribution (aka. inverse Beta)
#'
#' @source Bases on Bourguignon, M., Santos-Neto, M., & de Castro, M. (2018).
#' A new regression model for positive data (<https://arxiv.org/abs/1804.07734>)
#'
#' @details The beta prime distribution has density
#' \deqn{f(y) = \frac{y^{(\mu(\Phi+1)-1)} (1+y)^{(-(\mu(\Phi+1)+\Phi+2))}} {\Beta(\mu(1+\Phi), \Phi +2)}}
#' @details With the usual beta prime parameters
#' \deqn{\beta = \Phi + 2, \alpha = \mu(\Phi + 1)}
#'
#' @param x Value, x > 0.
#' @param mu Mean, mu > 0.
#' @param phi Precision, phi > 0.
#' @param log Optional argument. If true, returns the log density.
#'
#' @return Density of the pdf given x, mu and phi.
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 20, length.out = 1000)
#' plot(x, dbetaprime(x, mu = 4, phi = 2), type = "l")
dbetaprime <- function(x, mu, phi, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("beta prime is only defined for x > 0")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("beta prime is only defined for phi > 0")
  }
  if (isTRUE(any(mu <= 0))) {
    stop("beta prime is only defined for mu > 0")
  }

  # calculate the second argument for beta prime, given mu
  beta <- phi + 2
  alpha <- mu * (phi + 1)

  lpdf <- (alpha - 1) *
    log(x) +
    (-(alpha + beta)) * log1p(x) -
    lbeta(alpha, beta)

  # return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function of the beta prime distribution
#'
#' @param p Probabilities, for which to calculate the quantiles
#' @param mu Mean, mu > 0.
#' @param phi Precision, phi > 0.
#'
#' @return Quantiles of the beta prime distribution
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' plot(x, qbetaprime(x, mu = 1, phi = 2), type = "l")
qbetaprime <- function(p, mu, phi) {
  # check the arguments
  if (isTRUE(any(p < 0 | p > 1))) {
    stop("p has to be in an interval of [0, 1]")
  }
  if (isTRUE(phi <= 0)) {
    stop("beta prime is only defined for phi > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("beta prime is only defined for mu > 0")
  }

  # calculate argument alpha of phi/beta prime
  qb <- qbeta(p, mu * (phi + 1), phi + 2)

  # now calculate the qbetaprime using log rules log(qbeta) - log(1 - qbeta)
  lqbp <- log(qb) - log1p(-qb)
  return(exp(lqbp))
}

#' RNG for the beta prime distribution
#'
#' @param n Number of samples.
#' @param mu Mean, mu > 0.
#' @param phi Precision, phi > 0.
#'
#' @return Random numbers from the beta prime distribution.
#' @export
#'
#' @examples hist(rbetaprime(100, mu = 1, phi = 2))
rbetaprime <- function(n, mu = 1, phi = 1) {
  # check the arguments
  if (isTRUE(phi <= 0)) {
    stop("beta prime is only defined for phi > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("beta prime is only defined for mu > 0")
  }
  return(qbetaprime(runif(n, min = 0, max = 1), mu, phi))
}

#' Log-Likelihood of the beta prime distribution
#'
#' @param i brms indices
#' @param prep brms data
#'
#' @return Log-Likelihood of beta prime given data in prep
log_lik_betaprime <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  y <- prep$data$Y[i]
  return(dbetaprime(y, mu, phi, log = TRUE))
}


#' posterior_predict for the beta prime distribution
#'
#' @param i brms indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_betaprime <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  return(rbetaprime(prep$ndraws, mu, phi))
}

#' posterior_epred for the beta prime distribution
#'
#' @param prep brms data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_betaprime <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' Beta prime brms custom family
#'
#' @param link Link function for function
#' @param link_phi Link function for beta argument
#'
#' @return brms beta prime distribution family
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
betaprime <- function(link = "log", link_phi = "log") {
  family <- brms::custom_family(
    "betaprime",
    dpars = c("mu", "phi"),
    links = c(link, link_phi),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_betaprime,
    posterior_predict = posterior_predict_betaprime,
    posterior_epred = posterior_epred_betaprime
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real betaprime_lpdf(real y, real mu, real phi) {
        real beta = phi + 2;
        real alpha = mu * (phi + 1);
        return  (alpha-1) * log(y) +
                (-(alpha + beta)) * log1p(y) -
                lbeta(alpha, beta);
      }

      real betaprime_rng(real mu, real phi) {
        real rb = beta_rng(mu * (phi +1), phi + 2);
        return (rb / (1 - rb));
      }",
    block = "functions"
  )
  return(family)
}
