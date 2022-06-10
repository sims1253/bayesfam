#' Custom Beta distribution RNG
#'
#' @param n Number of draws.
#' @param mu Mean
#' @param phi Precision
#'
#' @return n samples beta distributed.
#' @export
#'
#' @examples samples <- rbeta_custom(1000, 0.5, 1)
#' hist(samples)
rbeta_custom <- function(n, mu, phi) {
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0,1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  return(rbeta(n, mu * phi, (1 - mu) * phi))
}

#' Custom Beta distribibution density.
#'
#' @param x x-value
#' @param mu Mean
#' @param phi Precision
#'
#' @details The beta-prime distribution has density
#' \deqn{f(y) = {y^{\mu \phi - 1}(1 - y)^{(1 - \mu) \phi - 1} } / beta(\mu \phi, (1 - \mu) \phi) }
#'
#' @return PDF of Custom Beta distribution, with mean parameterasation
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' phi <- 2
#' mean <- 0.5
#' y <- bayesim::dbeta_custom(x, mu = mean, phi = phi)
#' plot(x, y, type = "l", ylab = "Density", main = "dbeta_custom(mu = 0.5, phi = 2)")
#' # Compare to online ressources
dbeta_custom <- function(x, mu, phi) {
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0, 1).")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("P must be above 0.")
  }
  dbeta(x, mu * phi, (1 - mu) * phi)
}

#' Log-Likelihood vignette for the Custom-Beta distribution, with Mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of Beta-Custom given data in prep
#'
#' @examples
log_lik_beta <- function(i, prep) {
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  y <- prep$data$Y[i]
  return(dbeta_custom(y, mu, phi, log = TRUE))
}

#' Posterior predictionLog-Likelihood vignette for the Custom-Beta distribution, with Mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of Beta-Custom, given data in prep
#'
#' @examples
posterior_predict_beta <- function(i, prep, ...) {
  mu <- get_dpar(prep, "mu", i = i)
  phi <- get_dpar(prep, "phi", i = i)
  return(rbeta_custom(prep$ndraws, mu, phi))
}

#' Title
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
#'
#' @examples
posterior_epred_beta <- function(prep) {
  mu <- get_dpar(prep, "mu")
  return(mu)
}

#' Custom-Beta Stan-implementation in mean parametrization.
#'
#' @param link Link function for function
#' @param link_phi Link function for phi argument
#'
#' @return BRMS Beta-Custom distribution family
#'
#' @examples n <- 10000
#' a <- rnorm(n)
#' data <- list(a = a, y = bayesim::rbeta_custom(n, bayesim::inv_logit(0.5 * a + 1), 2))
#' fit1 <- brm(y ~ 1 + a,
#'   data = data, family = bayesim::beta_custom(),
#'   stanvars = bayesim::beta_custom()$stanvars, backend = "cmdstan"
#' )
#' plot(fit1)
beta_custom <- function(link = "logit", link_phi = "log") {
  family <- brms::custom_family(
    "beta_custom",
    dpars = c("mu", "phi"),
    links = c(link, link_phi),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_beta,
    posterior_predict = posterior_predict_beta,
    posterior_epred = posterior_epred_beta
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real beta_custom_lpdf(real y, real mu, real phi) {
        return(beta_lpdf(y | mu*phi, (1-mu)*phi));
      }

      real beta_custom_rng(real mu, real phi) {
        return(beta_rng(mu*phi, (1-mu)*phi));
      }",
    block = "functions"
  )
  return(family)
}
