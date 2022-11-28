#' Cauchitnormal Density distribution.
#'
#' @param x Value space of the distribution, x e (0, 1)
#' @param mu Median parameter, mu is already Cauchit-transformed, mu unbound
#' @param sigma Sigma shape parameter, sigma >= 0
#' @param log Bool argument, if true, returns the logarithmic density
#'
#' @return Normal Distribution Density with Cauchit link function
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, dcauchitnormal(x, mu = 0.5, sigma = 2), type = "l")
dcauchitnormal <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  logpdf <-
    -(log(sigma) + 0.5 * (log(2 * pi))) +
    log(2 * pi) - log(cos(2 * pi * (x - 0.5)) + 1) +
    -(cauchit(x) - mu)^2 / (2 * sigma^2)
  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Cauchitnormal RNG-function
#'
#' @param n Number of draws
#' @param mu Median parameter, mu unbound, mu already cauchit transformed
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @returns n chauchitnormal distributed samples
#'
#' @export
#'
#' @examples hist(rcauchitnormal(100, 0.5, 2))
rcauchitnormal <- function(n, mu, sigma) {
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    inv_cauchit(rnorm(n, mu, sigma))
  )
}

#' Log-Likelihood vignette for the Cauchitnormal Distribution
#'
#' @param i Indices
#' @param prep BRMS data
#'
#' @return log_likelihood of the Cauchitnormal Distribution, given some BRMS data.
log_lik_cauchitnormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dcauchitnormal(y, mu, sigma, log = TRUE))
}

#' Posterior-predict vignette for the Cauchitnormal Distribution
#'
#' @param i Indices
#' @param prep BRMS data
#' @param ...
#'
#' @return The posterior prediction of the Cauchitnormal Distribution, given some BRMS data.
posterior_predict_cauchitnormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rcauchitnormal(prep$ndraws, mu, sigma))
}

#' Posterior expected value prediction. Mean undefined for Cauchitnormal.
#'
#' @param prep BRMS data
#'
#' @return Nothing
posterior_epred_cauchitnormal <- function(prep) {
  # https://doi.org/10.1080/03610926.2020.1752723 might solve this
  stop("Due to the mean not having an analytical solution for the cauchit-normal
        distribution, posterior_epred is currently not supported.")
}

#' Custom BRMS family Cauchitnormal
#'
#' @param link Link function argument (as string) for Median argument. Left as identity!
#' @param link_sigma Link function argument (as string) for Shape argument
#'
#' @return Cauchitnormal BRMS model-object
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' a <- rnorm(1000)
#' data <- list(a = a, y = rcauchitnormal(1000, 0.5 * a + 1, 2))
#' # BBmisc::surpressAll necassary to keep the test output clean
#' BBmisc::suppressAll({
#'   fit1 <- brms::brm(y ~ 1 + a,
#'     data = data, family = cauchitnormal(),
#'     stanvars = cauchitnormal()$stanvars, backend = "cmdstanr", cores = 4
#'   )
#' })
#' plot(fit1)
cauchitnormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "cauchitnormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_cauchitnormal,
    posterior_predict = posterior_predict_cauchitnormal,
    posterior_epred = posterior_epred_cauchitnormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real cauchitnormal_lpdf(real y, real mu, real sigma) {
      return -(log(sigma) + 0.5 * (log(2 * pi()))) +
              log(2 * pi()) - log(cos(2 * pi() * (y - 0.5)) + 1) +
              -(tan(pi() * (y - 0.5)) - mu)^2 / (2 * sigma^2);
      }

      real cauchitnormal_rng(real mu, real sigma) {
        return cauchy_cdf(normal_rng(mu, sigma), 0, 1);
      }",
    block = "functions"
  )
  return(family)
}
