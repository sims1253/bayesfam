#' Probability density function for the Gompertz distribution, with Median parametrization.
#'
#' @param x Model space, defined for x >= 0
#' @param mu Median parameter of pdf, mu > 0
#' @param eta Second shape parameter of Gompertz, defined for eta > 1
#' @param log Optional argument. If TRUE, returns log(pdf). Normally False.
#'
#' @return PDF of gompertz distribution, with median parametrization.
#' @export
#'
#' @examples x <- seq(from = 0, to = 5, length.out = 100)
#' eta <- 0.1
#' b <- 1
#' median <- (1 / b) * log((-1 / eta) * log(1 / 2) + 1)
#' y <- bayesim::dgompertz(x, mu = median, eta = eta)
#' plot(x, y, type = "l", ylab = "Density", main = "dgompertz(mu=2.0708, eta=0.1) or dgompertz(b=1, eta=0.1)")
#' # Compare to online ressources
dgompertz <- function(x, mu, b, log = FALSE) {
  lpdf <- log(-b * log(0.5)) - log(exp(mu * b) - 1) + b * x - (-log(0.5)/(exp(mu * b) - 1)) * (exp(b * x) - 1)
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}


#' Quantile function for the Gompertz distribution, with Median parametrization.
#'
#' @param p Quantile to be calculated
#' @param mu Median argument of Gompertz
#' @param eta Eta argument of Gompertz
#'
#' @return Inverse of CDF, calculates a value, given a probability p
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' plot(x, bayesim::qgompertz(x, mu = 2, eta = 0.1), type = "l", ylab = "Quantile", main = "apex-after-origin Gompertz(mu=2, eta=0.1)")
qgompertz <- function(p, mu, b) {
  a <- -(b * log(0.5)) / (exp(mu * b) - 1)
  x <- (1 / b) * log1p(-(b/a) * log1p(-p))
  return(x)
}

#' RNG function for the Gompertz distribution, with Median parametrization.
#'
#' @param n Number of draws
#' @param mu Median argument of Gompertz
#' @param eta Eta argument of Gompertz
#'
#' @return A Gompertz distributed RNG vector of size n
#' @export
#'
#' @examples y <- bayesim::rgompertz(n, mu = 2, eta = 0.1)
#' hist(y, main = c(paste("Median:", median(y)), " for RNG of apex-after-origin Gompertz(mu=2, eta=0.1)"))
rgompertz <- function(n, mu, b) {
  return(qgompertz(runif(n), mu, b))
}

#' Log-Likelihood vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of gompertz given data in prep
#'
#' @examples
log_lik_gompertz <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  b <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  return(dgompertz(y, mu, b, log = TRUE))
}

#' Posterior-Prediction vignette for the Gompertz distribution, with Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of gompertz, given data in prep
#'
#' @examples
posterior_predict_gompertz <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  b <- brms::get_dpar(prep, "beta", i = i)
  return(rgompertz(prep$ndraws, mu, b))
}

#' Expectation-Predict vignette for the Gompertz distribution, with Median parametrization.
#' Not defined for the Gompertz family.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
#'
#' @examples
posterior_epred_gompertz <- function(prep) {
  stop("posterior_epred is not defined for the gompertz family")
}


#' Gompertz Stan-implementation in median parametrization.
#'
#' @param link Link function for function
#' @param link_eta Link function for eta argument
#'
#' @return BRMS gompertz distribution family
#' @export
#'
#' @examples data <- list(a = a, y = bayesim::rgompertz(n, exp(0.5 * a + 1), 0.2))
#' fit1 <- brm(y ~ 1 + a,
#'   data = data, family = bayesim::gompertz(),
#'   stanvars = bayesim::gompertz()$stanvars, backend = "cmdstan"
#' )
#' plot(fit1)
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
        real log_a = log(-beta * log(0.5)) - log(exp(mu * beta) - 1);
        real a_div_b = -log(0.5)/(exp(mu * beta) - 1);
        real lpdf = log_a + beta * y - a_div_b * (exp(beta * y) - 1);
        return(lpdf);
      }
      real gompertz_rng(real mu, real beta) {
        real a = -(beta * log(0.5)) / (exp(mu * beta) - 1);
        real x = (1 / beta) * log1p(-(beta/a) * log1p(uniform_rng(-1, 0)));
        return(x);
      }",
    block = "functions"
  )
  return(family)
}
