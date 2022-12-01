#' Logitnormal density distribution in median parametrization.
#'
#' @param x Value space of the distribution, x e (0, 1)
#' @param mu Median parameter, mu is already logit-transformed, mu unbound
#' @param sigma Sigma shape parameter, sigma >= 0
#' @param log Bool argument, if true, returns the logarithmic density
#'
#' @return Normal distribution density with logit link function
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, dlogitnormal(x, mu = 0.5, sigma = 2), type = "l")
dlogitnormal <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  logpdf <-
    -(log(sigma) + 0.5 * (log(2 * pi))) +
    (-(log(x) + log1p(-x))) +
    (-(logit(x) - mu)^2 / (2 * sigma^2))
  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Logitnormal RNG-function in median parametrization.
#'
#' @param n Number of draws
#' @param mu Median parameter, mu unbound, mu already logit transformed
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @returns n Logitnormal distributed samples
#'
#' @export
#'
#' @examples hist(rlogitnormal(100, 0.5, 2))
rlogitnormal <- function(n, mu, sigma) {
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    inv_logit(rnorm(n, mu, sigma))
  )
}

#' Log-Likelihood vignette for the Logitnormal distribution, in Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#'
#' @return log_likelihood of the Logitnormal distribution, given some BRMS data.
log_lik_logitnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dlogitnormal(y, mu, sigma, log = TRUE))
}

#' Posterior-predict vignette for the Logitnormal distribution, with Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#' @param ...
#'
#' @return The posterior prediction of the Logitnormal distribution, given some BRMS data.
posterior_predict_logitnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rlogitnormal(prep$ndraws, mu, sigma))
}

#' Posterior expected value prediction vignette for Logitnormal distribution.
#'
#' @param prep BRMS data
#'
#' @return Mean of Posterior
posterior_epred_logitnormal <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  return(exp(mu + sigma^2 / 2))
}

#' Custom BRMS family Logit-Normal in median parametrization.
#'
#' @param link Link function argument (as string) for Median argument. Left as identity!
#' @param link_sigma Link function argument (as string) for Shape argument
#'
#' @return Logitnormal BRMS model-object
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' a <- rnorm(1000)
#' data <- list(a = a, y = rlogitnormal(1000, 0.5 * a + 1, 2))
#' # refresh = 0 supresses chain updates
#' fit1 <- brms::brm(y ~ 1 + a, data = data,
#'  family = logitnormal(), stanvars = logitnormal()$stanvars,
#'  refresh = 0)
#' plot(fit1)
logitnormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "logitnormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_logitnormal,
    posterior_predict = posterior_predict_logitnormal,
    posterior_epred = posterior_epred_logitnormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real logitnormal_lpdf(real y, real mu, real sigma) {
        return -(log(sigma) + 0.5 * (log(2*pi()))) +
              -(log(y) + log1m(y)) +
              -((log(y) - log1m(y)) - mu)^2 / (2 * sigma^2);
      }

      real logitnormal_rng(real mu, real sigma) {
        return inv_logit(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
