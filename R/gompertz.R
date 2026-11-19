#' Probability density function for the Gompertz distribution, with Median parametrization.
#'
#' @param x Value
#' @param mu Median parameter
#' @param beta Scale parameter
#' @param log Optional argument. If TRUE, returns log(pdf).
#'
#' @details PDF of Gompertz implementation, with constant b:
#' \deqn{b(\mu,\eta) := (1 / \mu) * log1p((-1 / \eta) * log(0.5))}
#' \deqn{f(x) = \eta*b*exp(\eta + bx - \eta * e^{bx})}
#'
#' @return f(x | mu, eta)
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 5, length.out = 100)
#' plot(x, dgompertz(x, mu = 2, beta = 4), type = "l")
dgompertz <- function(x, mu, beta, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("dgompertz is only defined for x > 0")
  }
  if (isTRUE(any(mu <= 0))) {
    stop("dgompertz is only defined for mu > 0")
  }
  if (isTRUE(any(beta <= 0))) {
    stop("dgompertz is only defined for beta > 0")
  }

  lpdf <- log(-beta * log(0.5)) -
    log(expm1(mu * beta)) +
    (beta * x + (log(0.5) / expm1(mu * beta)) * expm1(beta * x))

  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function for the Gompertz distribution, with Median parametrization.
#'
#' @param p Quantile to be calculated
#' @param mu Median parameter
#' @param beta Scale parameter
#'
#' @return Inverse of CDF, calculates a value, given a probability p
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 100)
#' plot(x, qgompertz(x, mu = 10, beta = 1), type = "l")
qgompertz <- function(p, mu, beta) {
  # check the arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("qgompertz is only defined for 0 < p < 1")
  }
  if (isTRUE(any(mu <= 0))) {
    stop("qgompertz is only defined for mu > 0")
  }
  if (isTRUE(any(beta <= 0))) {
    stop("qgompertz is only defined for beta > 0")
  }
  return(
    log1p((expm1(mu * beta) / (log(0.5))) * log1p(-p)) / beta
  )
}

#' RNG function for the Gompertz distribution, with Median parametrization.
#'
#' @param n Number of draws
#' @param mu Median parameter
#' @param beta Scale parameter
#'
#' @return A Gompertz distributed RNG vector of size n
#' @export
#'
#' @examples hist(rgompertz(n = 100, mu = 2, beta = 0.1))
rgompertz <- function(n, mu = 1, beta = 0.5) {
  if (isTRUE(mu <= 0)) {
    stop("rgompertz is only defined for mu > 0")
  }
  if (isTRUE(beta <= 0)) {
    stop("rgompertz is only defined for beta > 0")
  }
  return(
    log1p((expm1(mu * beta) / (log(0.5))) * log(runif(n))) / beta
  )
}


#' Log-Likelihood vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i brms indices
#' @param prep brms data
#'
#' @return Log-Likelihood of gompertz given data in prep
log_lik_gompertz <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dgompertz(y, mu, beta, log = TRUE))
}

#' Posterior-Prediction vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i brms indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Posterior prediction of gompertz, given data in prep
posterior_predict_gompertz <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  return(rgompertz(prep$ndraws, mu, beta))
}

#' Scaled exponential integral exp(x) * E1(x) for x > 0.
#'
#' Series for x <= 1, continued fraction for x > 1; relative error about
#' 1e-12.
#'
#' @param x Positive values
#'
#' @return exp(x) * E1(x), same dimensions as `x`.
#' @noRd
exp_E1 <- function(x) {
  out <- x
  small <- x <= 1
  # E1(x) = -gamma - log(x) - sum_k (-x)^k / (k * k!)
  xs <- x[small]
  total <- 0
  term <- 1
  for (k in 1:30) {
    term <- -term * xs / k
    total <- total + term / k
  }
  out[small] <- exp(xs) * (-0.5772156649015329 - log(xs) - total)
  # exp(x) E1(x) = 1 / (x + 1 - 1 / (x + 3 - 4 / (x + 5 - ...)))
  xl <- x[!small]
  f <- xl + 201
  for (k in 100:1) {
    f <- xl + 2 * k - 1 - k^2 / f
  }
  out[!small] <- 1 / f
  return(out)
}

#' Mean of the Gompertz distribution, with Median parametrization.
#'
#' With eta = log(2) / expm1(mu * beta), the survival function is
#' S(x) = exp(-eta * (exp(beta * x) - 1)), so the mean is
#' integral S(x) dx = exp(eta) * E1(eta) / beta.
#'
#' @param mu Median parameter
#' @param beta Scale parameter
#'
#' @return Mean, same dimensions as `mu`.
#' @noRd
gompertz_mean <- function(mu, beta) {
  x <- mu * beta
  out <- exp_E1(log(2) / expm1(x))
  # where z = log(2) / expm1(x) underflows, exp(z) E1(z) is
  # -gamma - log(z) + O(z log(z)), with log(z) finite even if expm1(x) is not
  lz <- log(log(2)) - x - log(-expm1(-x))
  tiny <- lz < log(.Machine$double.xmin)
  out[tiny] <- -0.5772156649015329 - lz[tiny]
  return(out / beta)
}

#' Expectation-Predict vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param prep brms data
#'
#' @return Mean of the posterior predictive distribution, draws x observations
posterior_epred_gompertz <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  beta <- brms::get_dpar(prep, "beta")
  return(gompertz_mean(mu, beta))
}


#' Custom Gompertz brms-implementation in median parametrization.
#'
#' @param link Link function for function
#' @param link_b Link function for eta argument
#'
#' @return brms gompertz distribution family
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(a = a, y = rgompertz(1000, mu = exp(0.5 * a + 1), beta = 0.1))
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = gompertz(), stanvars = gompertz()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
gompertz <- function(link = "log", link_b = "log") {
  family <- brms::custom_family(
    "gompertz",
    dpars = c("mu", "beta"),
    links = c(link, link_b),
    lb = c(0, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_gompertz,
    posterior_predict = posterior_predict_gompertz,
    posterior_epred = posterior_epred_gompertz
  )

  family$stanvars <- brms::stanvar(
    scode = "
      real gompertz_lpdf(real y, real mu, real beta) {
        return(
          log(-beta * log(0.5)) -
          log(expm1(mu * beta)) +
          (beta * y + (log(0.5) / expm1(mu * beta)) * expm1(beta * y))
        );
      }
      real gompertz_rng(real mu, real beta) {
      return(
        log1p((expm1(mu * beta) / (log(0.5))) * log(uniform_rng(0,1))) / beta
      );
      }",
    block = "functions"
  )
  return(family)
}
