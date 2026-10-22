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
#' @param n number of observations
#' @param mu Median parameter, mu unbound, mu already logit transformed
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @returns n Logitnormal distributed samples
#'
#' @export
#'
#' @examples hist(rlogitnormal(100, 0.5, 2))
rlogitnormal <- function(n, mu = 0, sigma = 1) {
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  return(
    inv_logit(rnorm(n, mu, sigma))
  )
}

#' Log-Likelihood vignette for the Logitnormal distribution, in Median parametrization.
#'
#' @param i Indices
#' @param prep brms data
#'
#' @return log_likelihood of the Logitnormal distribution, given some brms data.
log_lik_logitnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dlogitnormal(y, mu, sigma, log = TRUE))
}

#' Posterior-predict vignette for the Logitnormal distribution, with Median parametrization.
#'
#' @param i Indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return The posterior prediction of the Logitnormal distribution, given some brms data.
posterior_predict_logitnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rlogitnormal(prep$ndraws, mu, sigma))
}

#' Expected response of the logit-normal distribution.
#'
#' Computes E[plogis(Z)] for Z ~ N(mu, sigma) (the expectation that
#' `posterior_epred` promises) with trapezoid sums on fixed nodes, which are
#' accurate to an absolute error of about 1e-14 for all mu and sigma.
#'
#' @param mu Median parameter, unbound
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @return Expected response E[plogis(Z)], same dimensions as `mu`.
#' @noRd
logitnormal_mean <- function(mu, sigma) {
  sigma <- sigma + 0 * mu
  mu <- mu + 0 * sigma
  out <- 0 * mu
  # sigma <= 1: integrate plogis(mu + sigma * t) over t ~ N(0, 1)
  narrow <- sigma <= 1
  m <- mu[narrow]
  s <- sigma[narrow]
  total <- 0
  for (t in seq(-9, 9, by = 0.5)) {
    total <- total + 0.5 * dnorm(t) * plogis(m + s * t)
  }
  out[narrow] <- total
  # sigma > 1 makes that integrand too steep, so use the dual form: with
  # L ~ Logistic(0, 1), E[plogis(Z)] = P(L < Z) = E[pnorm((mu + L) / sigma)]
  m <- mu[!narrow]
  s <- sigma[!narrow]
  total <- 0
  for (l in seq(-36, 36, by = 0.5)) {
    total <- total + 0.5 * dlogis(l) * pnorm((m + l) / s)
  }
  out[!narrow] <- total
  return(out)
}

#' Posterior expected value prediction vignette for Logitnormal distribution.
#'
#' @param prep brms data
#'
#' @return Mean of Posterior
posterior_epred_logitnormal <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  return(logitnormal_mean(mu, sigma))
}

#' Custom brms family Logit-Normal in median parametrization.
#'
#' @param link Link function argument (as string) for Median argument. Left as identity!
#' @param link_sigma Link function argument (as string) for Shape argument
#'
#' @return Logitnormal brms model-object
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(a = a, y = rlogitnormal(1000, 0.5 * a + 1, 2))
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = logitnormal(), stanvars = logitnormal()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
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
