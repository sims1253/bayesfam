#' symlognormal density function
#'
#' @source Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
#'         Mastering Diverse Domains through World Models.
#'         (<https://doi.org/10.48550/arXiv.2301.04104>)
#'
#' @param x Value space of the distribution, x unbound
#' @param mu Median parameter, mu unbound
#' @param sigma Sigma shape parameter, sigma >= 0
#' @param log Bool argument, if true, returns the logarithmic density
#'
#' @return Normal distribution density with symlog link function
#' @export
#'
#' @examples x <- seq(from = -5, to = 5, length.out = 1000)
#' plot(x, dsymlognormal(x, mu = -0.2, sigma = 0.4), type = "l")
dsymlognormal <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }

  symlog_jacobian_adjustment <- ifelse(
    x == 0,
    0,
    log(x * sign(x)) - log((abs(x) + x^2))
  )
  logpdf <-
    -(log(sigma) + 0.5 * (log(2 * pi))) +
    symlog_jacobian_adjustment +
    (-(symlog(x) - mu)^2 / (2 * sigma^2))

  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' symlognormal RNG-function
#'
#' @source Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
#'         Mastering Diverse Domains through World Models.
#'         (<https://doi.org/10.48550/arXiv.2301.04104>)
#'
#' @param n Number of draws
#' @param mu Median parameter, mu unbound, mu already symlog transformed
#' @param sigma Sigma shape parameter, sigma > 0
#'
#' @returns n symlognormal distributed samples
#'
#' @export
#'
#' @examples hist(rsymlognormal(100, 0.5, 2))
rsymlognormal <- function(n, mu = 0, sigma = 1) {
  if (isTRUE(any(sigma < 0))) {
    stop("sigma must be above or equal to 0.")
  }
  return(
    inv_symlog(rnorm(n, mu, sigma))
  )
}

#' Log-Likelihood for the symlognormal distribution,
#'
#' @param i Indices
#' @param prep brms data
#'
#' @return log_likelihood of the symlognormal distribution, given some brms data.
log_lik_symlognormal <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dsymlognormal(y, mu, sigma, log = TRUE))
}

#' Posterior-predict for the symlognormal distribution, with Median parametrization.
#'
#' @param i Indices
#' @param prep brms data
#' @param ... catchall
#'
#' @return The posterior prediction of the symlognormal distribution, given some brms data.
posterior_predict_symlognormal <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rsymlognormal(prep$ndraws, mu, sigma))
}

#' Posterior epred for symlognormal distribution.
#'
#' @param prep brms data
#'
#' @return warning
posterior_epred_symlognormal <- function(prep) {
  warning(
    "posterior_epred is not defined for the symlog normal as I don't know a mean formula."
  )
}

#' Custom brms family symlog-Normal
#'
#' @source Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
#'         Mastering Diverse Domains through World Models.
#'         (<https://doi.org/10.48550/arXiv.2301.04104>)
#'
#' @param link Link function argument (as string) for Median argument. Left as identity!
#' @param link_sigma Link function argument (as string) for Shape argument
#'
#' @return symlognormal brms model-object
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(a = a, y = rsymlognormal(1000, 0.5 * a + 1, 2))
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = symlognormal(), stanvars = symlognormal()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
symlognormal <- function(link = "identity", link_sigma = "log") {
  stopifnot(link == "identity")
  family <- brms::custom_family(
    "symlognormal",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(-NA, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_symlognormal,
    posterior_predict = posterior_predict_symlognormal,
    posterior_epred = posterior_epred_symlognormal
  )
  family$stanvars <- stanvars <- brms::stanvar(
    scode = "
      int sign(real x) {
        if (x > 0) {
          return 1;
        } else if (x < 0) {
          return -1;
        } else {
          return 0;
        }
      }

      real symlog(real x) {
        return(sign(x) * log1p(abs(x)));
      }

      real inv_symlog(real x) {
        return sign(x)*(expm1(abs(x)));
      }

      real symlognormal_lpdf(real y, real mu, real sigma) {
        real symlog_jacobian_adjustment = 0;
        if (y == 0) {
          symlog_jacobian_adjustment = 0;
        } else {
          symlog_jacobian_adjustment = log(y * sign(y)) - log((abs(y)+y^2));
        }

        return -(log(sigma) + 0.5 * (log(2*pi()))) +
              symlog_jacobian_adjustment +
              (-(symlog(y) - mu)^2 / (2 * sigma^2));
      }

      real symlognormal_rng(real mu, real sigma) {
        return inv_symlog(normal_rng(mu, sigma));
      }",
    block = "functions"
  )
  return(family)
}
