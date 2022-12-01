#' Kumaraswamy density function in median parametrisation.
#'
#' @param x value space, x e (0, 1)
#' @param mu Median parameter of pdf, mu e (0, 1)
#' @param p shape parameter, p > 0
#' @param log if true, returns log(pdf). Normally FALSE.
#'
#' @details \deqn{q(\mu, p) = -\frac{log(2)}{log(1-\mu^p)}}
#' @details \deqn{f(y | \mu, p) = pqx^{p-1}(1-x^p)^{q-1}}
#'
#' @return f(x | mu, p)
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, dkumaraswamy(x, mu = 0.5, p = 2), type = "l")
dkumaraswamy <- function(x, mu, p, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop("The mean must be in (0,1).")
  }
  mu[which(mu > 0.999999)] <- 0.999999
  mu[which(mu < 0.000001)] <- 0.000001
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  logpdf <- log(p) +
    log(log(2)) -
    log(-(log1p(-mu^p))) +
    (p - 1) * log(x) +
    ((-(log(2) / log1p(-mu^p))) - 1) * log1p(-x^p)

  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Kumaraswamy RNG function in Median parametrization.
#'
#' @param n Number of samples to draw, as a natural number scalar.
#' @param mu Median parameter, mu e (0, 1)
#' @param p Phi shape parameter, Phi > 0
#'
#' @return n samples in Kumaraswamy distribution.
#' @export
#'
#' @examples hist(rkumaraswamy(10000, mu = 0.5, p = 4))
rkumaraswamy <- function(n, mu, p) {
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mean must be in (0,1).")
  }
  mu[which(mu > 0.999999)] <- 0.999999
  mu[which(mu < 0.000001)] <- 0.000001
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  q <- -(log(2) / log1p(-mu^p))
  return(
    (1 - (1 - runif(n, min = 0, max = 1))^(1 / q))^(1 / p)
  )
  # TODO: Kumaraswamy BRMS does not like RNG values at the boundary.
  # maybe one might also clip them in the RNG return already?
  # (Instead of clipping them before feeding the BRMS)
}

#' Quantile function of the Kumaraswamy distribution in Median parametrisation.
#'
#' @param u Quantile to be calculated, u e (0, 1)
#' @param mu Median parameter, mu e (0, 1)
#' @param p Phi shape parameter, p > 0
#'
#' @return q(u | mu, p)
#' @export
#'
#' @examples u <- seq(from = 0.01, to = 0.09, length.out = 1000)
#' plot(u, qkumaraswamy(u, mu = 0.5, p = 2), type = "l")
qkumaraswamy <- function(u, mu = 0.5, p = 1) {
  if (isTRUE(any(u <= 0 | u >= 1))) {
    stop("u must be in (0,1).")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  return(
    (1 -
      (1 - u)^(1 /
        (-(log(2) / log1p(-mu^p))))
    )^(1 / p)
  )
}

#' Kumaraswamy CDF in median parametrisation
#'
#' @param x CDF of x over lower tail, x e (0, 1)
#' @param mu Median parameter, mu e (0, 1)
#' @param p shape parameter, p > 0
#'
#' @return p(x | mu, p)
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, pkumaraswamy(x, mu = 0.5, p = 1), type = "l")
pkumaraswamy <- function(x, mu = 0.5, p = 1) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  q <- -(log(2) / log1p(-mu^p))
  return(1 + (x^p - 1)^q)
}

#' Log-Likelihood vignette for the Kumaraswamy distribution, in Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of Kumaraswamy given data in prep
log_lik_kumaraswamy <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  y <- prep$data$Y[i]
  return(dkumaraswamy(y, mu, p, log = TRUE))
}

#' Posterior prediction vignette for the Kumaraswamy distribution, in Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of Kumaraswamy, given data in prep
posterior_predict_kumaraswamy <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  return(rkumaraswamy(prep$ndraws, mu, p))
}

#' Posterior expected value prediction of the Kumaraswamy implementation.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
posterior_epred_kumaraswamy <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  p <- brms::get_dpar(prep, "p")
  q <- -(log(2) / log1p(-mu^p))
  return(q * beta((1 + 1 / p), q))
}

#' Kumaraswamy BRMS-implementation in median parametrization.
#'
#' @param link Link function for function
#' @param link_p Link function for p argument
#'
#' @return BRMS Beta-Custom distribution family
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' a <- rnorm(1000)
#' data <- list(a = a, y = rkumaraswamy(1000, brms::inv_logit_scaled(0.5 * a + 1), 2))
#' # BBmisc::surpressAll necassary to keep the test output clean
#' BBmisc::suppressAll({
#'   fit1 <- brms::brm(y ~ 1 + a,
#'     data = data, family = kumaraswamy(),
#'     stanvars = kumaraswamy()$stanvars
#'   )
#' })
#' plot(fit1)
kumaraswamy <- function(link = "logit", link_p = "log") {
  family <- brms::custom_family(
    "kumaraswamy",
    dpars = c("mu", "p"),
    links = c(link, link_p),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_kumaraswamy,
    posterior_predict = posterior_predict_kumaraswamy,
    posterior_epred = posterior_epred_kumaraswamy
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real kumaraswamy_lpdf(real y, real mu, real p) {
         return  (log(p) + log(log(2)) - log(-(log1m(mu^p))) + (p-1) * log(y) +
                 ((-(log(2)/log1m(mu^p)))-1) * log1m(y^p));
      }

      real kumaraswamy_rng(real mu, real p) {
         return ((1-(1-uniform_rng(0, 1))^(1/(-(log(2)/log1m(mu^p)))))^(1/p));
      }",
    block = "functions"
  )
  return(family)
}
