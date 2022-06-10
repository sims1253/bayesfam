#' Probability density function for the Beta-Prime distribution (aka. inverse Beta)
#'
#' @source Bases on Bourguignon, M., Santos-Neto, M., & de Castro, M. (2018).
#' A new regression model for positive data (https://arxiv.org/abs/1804.07734)
#'
#' @details The beta-prime distribution has density
#' \deqn{f(y) = y^{(\mu(\Phi+1)-1)} (1+y)^{(-(\mu(\Phi+1)+\Phi+2))} / beta(\mu(1+\Phi), \Phi +2)}
#'
#' @param x Value, x > 0.
#' @param mu Mean, mu > 0.
#' @param phi Precision, phi > 0.
#' @param log Optional argument. If true, returns the log density.
#'
#' @return density of the pdf given x, mu and phi.
#' @export
#'
#' @examples x <- seq(from = 0, to = 100, length.out = 1000)
#' y <- dbetaprime(x, mu = 4, phi = 2)
#' plot(x, y, type = "l", ylab = "Density", main = "dbetaprime(mu=4, phi=2)")
dbetaprime <- function(x, mu, phi, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("betaprime is only defined for x > 0")
  }
  if (isTRUE(phi <= 0)) {
    stop("betaprime is only defined for phi > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("betaprime is only defined for mu > 0")
  }

  # calculate the second argument for beta-prime, given mu
  beta <- phi + 2
  alpha <- mu * (phi + 1)

  lpdf <- (alpha - 1) * log(x) +
    (-(alpha + beta)) * log1p(x) -
    log(beta(alpha, beta))

  # return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Title
#'
#' @param p Probabilities, for which to calculate the quantiles
#' @param mu Mean, mu > 0.
#' @param phi Precision, phi > 0.
#'
#' @return Quantiles of the beta-prime distribution, given p, mu and phi
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#  y = bayesim::qbetaprime(x, mu = 1, phi = 2)
#  plot(x, y, type="l", ylab = "Quantile", main = "left-leaning Beta-Prime(mu=1,phi=2)"))
qbetaprime <- function(p, mu, phi) {
  # check the arguments
  if (isTRUE(any(p < 0 | p > 1))) {
    stop("p has to be in an interval of [0, 1]")
  }
  if (isTRUE(phi <= 0)) {
    stop("betaprime is only defined for phi > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("betaprime is only defined for mu > 0")
  }

  # calculate argument alpha of phi/betaprime
  qb <- qbeta(p, mu * (phi + 1), phi + 2)

  # now calculate the qbetaprime using log rules log(qbeta) - log(1 - qbeta)
  lqbp <- log(qb) - log1p(-qb)
  return(exp(lqbp))
}

#' Title
#'
#' @param n Number of beta-prime samples.
#' @param mu Mean, mu > 0.
#' @param phi Precision, phi > 0.
#'
#' @return Random numbers from the beta-prime distribution.
#' @export
#'
#' @examples y <- bayesim::rbetaprime(100, mu = 1, phi = 2)
#  hist(y, main = c(paste("Mean:", mean(y)), " for RNG of left-leaning Beta-Prime(mu=1,phi=2)"))
rbetaprime <- function(n, mu, phi) {
  # check the arguments
  if (isTRUE(phi <= 0)) {
    stop("betaprime is only defined for phi > 0")
  }
  if (isTRUE(mu <= 0)) {
    stop("betaprime is only defined for mu > 0")
  }
  return(qbetaprime(runif(n, min = 0, max = 1), mu, phi))
}

#' Title
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of betaprime given data in prep
#'
#' @examples
log_lik_betaprime <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  y <- prep$data$Y[i]
  return(dbetaprime(y, mu, phi, log = TRUE))
}


#' Title
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of beta-prime, given data in prep
#'
#' @examples
posterior_predict_betaprime <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  return(rbetaprime(prep$ndraws, mu, phi))
}

#' Title
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
#'
#' @examples
posterior_epred_betaprime <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' Title
#'
#' @param link Link function for function
#' @param link_beta Link function for beta argument
#'
#' @return BRMS beta-prime distribution family
#' @export
#'
#' @examples data <- list(a = a, y = bayesim::rbetaprime(
#'   n,
#'   exp(0.5 * rnorm(1000) + 1), 2
#' ))
#' fit1 <- brm(y ~ 1 + a,
#'   data = data, family = bayesim::betaprime(),
#'   stanvars = bayesim::betaprime()$stanvars, backend = "cmdstan"
#' )
#' plot(fit1)
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
                log(beta(alpha, beta));
      }

      real betaprime_rng(real mu, real phi) {
        real rb = beta_rng(mu * (phi +1), phi + 2);
        return (rb / (1 - rb));
      }",
    block = "functions"
  )
  return(family)
}
