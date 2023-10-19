#' Density of the Unit Lindley likelihood
#'
#' @param x vector of values, x e (0, 1)
#' @param mu Mean, mu e (0, 1)
#' @param log logical; if TRUE, log(pdf) is returned
#'
#' @return f(x | mu)
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 1000)
#' plot(x, dunit_lindley(x, 0.5))
dunit_lindley <- function(x, mu, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("The x argument has to be in the interval (0, 1)")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mu argument has to be in the interval (0, 1)")
  }

  lpdf <- 2 * log1p(-mu) - log(mu) - 3 * log1p(-x) - x * (1 - mu) / (mu * (1 - x))
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function of Unit-Lindley distribution
#'
#' @param p vector of probabilities
#' @param mu Mean, mu e (0, 1)
#'
#' @return q(p | mu)
#' @export
#'
#' @examples p <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(p, qunit_lindley(p, mu = 0.5))
qunit_lindley <- function(p, mu) {
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("The p argument has to be in the interval (0, 1)")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mu argument has to be in the interval (0, 1)")
  }

  lambert_term <- lamW::lambertWm1((1 / mu) * (p - 1) * exp(-1 / mu))
  return((1 / mu + lambert_term) / (1 + lambert_term))
}

#' Unit-Lindley RNG function
#'
#' @param n number of observations
#' @param mu Mean, mu e (0, 1)
#'
#' @return n samples drawn from the Unit Lindley distribution
#' @export
#'
#' @examples hist(runit_lindley(100, 0.5))
runit_lindley <- function(n, mu) {
  return(qunit_lindley(runif(n), mu))
}

#' Log-likelihood of the Unit Lindley family
#'
#' @param i Indices
#' @param prep brms data
#'
#' @return Log-Likelihood of brms data
log_lik_unitlindley <- function(i, prep) {
  return(
    dunit_lindley(
      x = prep$data$Y[i],
      mu = brms::get_dpar(prep, "mu", i = i),
      log = TRUE
    )
  )
}

#' Posterior predictions of the Unit Lindley family
#'
#' @param i Indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Draws from the posterior predictive distribution
posterior_predict_unitlindley <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  return(runit_lindley(prep$ndraws, mu))
}

#' Expectations of posterior predictions of the Unit Lindley family
#'
#' @param prep brms data
#'
#' @return Mean of the posterior predictive distribution
posterior_epred_unitlindley <- function(prep) {
  return(prep$dpars$mu)
}

#' Unit Lindley brms family
#'
#' @param link link for mu, default = logit
#'
#' @return Unit Lindley brms family
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = runit_lindley(n = 1000, mu = inv_logit(0.5 * a + 1)))
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = unit_lindley(), stanvars = unit_lindley()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
unit_lindley <- function(link = "logit") {
  family <- brms::custom_family(
    "unit_lindley",
    dpars = c("mu"),
    links = c(link),
    lb = c(0),
    ub = c(1),
    type = "real",
    log_lik = log_lik_unitlindley,
    posterior_predict = posterior_predict_unitlindley,
    posterior_epred = posterior_epred_unitlindley
  )

  family$stanvars <- brms::stanvar(
    scode = "
      real unit_lindley_lpdf(real y, real mu) {
        return 2*log1m(mu)-log(mu)-3*log1m(y)-y*(1-mu)/(mu*(1-y));
      }
    ",
    block = "functions"
  )
  return(family)
}
