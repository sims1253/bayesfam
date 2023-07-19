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
  if (isTRUE(mu <= 0)) {
    stop("dgompertz is only defined for mu > 0")
  }
  if (isTRUE(beta <= 0)) {
    stop("dgompertz is only defined for beta > 0")
  }

  lpdf <- log(-beta * log(0.5)) -
    log(exp(mu * beta) - 1) +
    (beta * x + (log(0.5) / (exp(mu * beta) - 1)) * (exp(beta * x) - 1))

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
  if (isTRUE(mu <= 0)) {
    stop("qgompertz is only defined for mu > 0")
  }
  if (isTRUE(beta <= 0)) {
    stop("qgompertz is only defined for beta > 0")
  }
  return(
    log1p(((exp(mu * beta) - 1) / (log(0.5))) * log1p(-p)) / beta
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
    log1p(((exp(mu * beta) - 1) / (log(0.5))) * log(runif(n))) / beta
  )
}


#' Log-Likelihood vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
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
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ... catchall
#'
#' @return Posterior prediction of gompertz, given data in prep
posterior_predict_gompertz <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  return(rgompertz(prep$ndraws, mu, beta))
}

#' Expectation-Predict vignette for the Gompertz distribution, with Median parametrization.
#' Not defined for the Gompertz family.
#'
#' @param prep BRMS data
#'
#' @return Nothing
posterior_epred_gompertz <- function(prep) {
  stop("posterior_epred is not defined for the gompertz family")
}


#' Custom Gompertz BRMS-implementation in median parametrization.
#'
#' @param link Link function for function
#' @param link_b Link function for eta argument
#'
#' @return BRMS gompertz distribution family
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(a = a, y = rgompertz(1000, mu = exp(0.5 * a + 1), beta = 0.1))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'  family = gompertz(), stanvars = gompertz()$stanvars,
#'  refresh = 0)
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
          log(exp(mu * beta) - 1) +
          (beta * y + (log(0.5) / (exp(mu * beta) - 1)) * (exp(beta * y) - 1))
        );
      }
      real gompertz_rng(real mu, real beta) {
      return(
        log1p(((exp(mu * beta) - 1) / (log(0.5))) * log(uniform_rng(0,1))) / beta
      );
      }",
    block = "functions"
  )
  return(family)
}
