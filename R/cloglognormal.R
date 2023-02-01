#' Probability Density function of the Cloglognormal-Distribution
#'
#' @param x Value space of the function, x e (0, 1)
#' @param mu Median parameter, mu is already Cloglog-transformed, mu unbound
#' @param sigma Shape parameter, sigma >= 0
#' @param log optional argument. If true, returns logarathmic probability. Default = FALSE
#'
#' @return Normal Distribution density with Cloglog link function
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, dcloglognormal(x, mu = 0.5, sigma = 2), type = "l")
dcloglognormal <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  logpdf <-
    -(log(sigma) + 0.5 * (log(2 * pi))) +
    -(log((x - 1) * log1p(-x))) +
    -(cloglog(x) - mu)^2 / (2 * sigma^2)
  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Cloglognormal RNG-function
#'
#' @param n Number of draws
#' @param mu Median parameter, mu unbound, mu already cloglog transformed
#' @param sigma Shape parameter
#'
#' @return n cloglog-normally distributed samples
#' @export
#'
#' @examples hist(rcloglognormal(100, 0.5, 2))
rcloglognormal <- function(n, mu = -0.36, sigma = 0.75) {
  if (isTRUE(any(sigma < 0))) {
    stop("P must be above or equal to 0.")
  }
  return(
    inv_cloglog(rnorm(n, mu, sigma))
  )
}

#' Log-Likelihood vignette for the Chauchitnormal distribution, with Median parametrization.
#'
#' @param i indices
#' @param prep BRMS data
#'
#' @return log_lik
log_lik_cloglognormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dcloglognormal(y, mu, sigma, log = TRUE))
}

#' Posterior-predict vignette for the Chauchitnormal distribution, with Median parametrization.
#'
#' @param i Indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of the data
posterior_predict_cloglognormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rcloglognormal(prep$ndraws, mu, sigma))
}

#' Posterior expected value prediction. Mean undefined for CLogLog-Normal
#'
#' @param prep BRMS data
#'
#' @return Nothing
posterior_epred_cloglognormal <- function(prep) {
  # https://doi.org/10.1080/03610926.2020.1752723 might solve this
  stop("Due to the mean not having an analytical solution for the cloglog-normal
        distribution, posterior_epred is currently not supported.")
}

#' Custom BRMS family CLogLog-Normal in median parametrization.
#'
#' @param link Link function argument (as string) for Median argument. Left as identity!
#' @param link_sigma Link function argument (as string) for Shape argument
#'
#' @return Cloglog BRMS model-object
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' cloglog_data <- rcloglognormal(1000, 0.5, 2)
#' # cloglognormal does not like values to close to the boundary
#' cloglog_data <- limit_data(cloglog_data, c(1e-12, 1 - 1e-12))
#' # BBmisc::surpressAll necassary to keep the test output clean
#' BBmisc::suppressAll({
#'   fit1 <- brms::brm(y ~ 1,
#'     data = list(y = cloglog_data), family = cloglognormal(),
#'     stanvars = cloglognormal()$stanvars, backend = "cmdstanr", cores = 4
#'   )
#' })
#' plot(fit1)
cloglognormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "cloglognormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_cloglognormal,
    posterior_predict = posterior_predict_cloglognormal,
    posterior_epred = posterior_epred_cloglognormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      real cloglognormal_lpdf(real y, real mu, real sigma) {
        return  -(log(sigma) + 0.5 * (log(2 * pi()))) +
                -(log((y - 1) * log1m(y))) +
                -(log(-log(1 - y)) - mu)^2 / (2 * sigma^2);
      }

      real cloglognormal_rng(real mu, real sigma) {
        return inv_cloglog(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
