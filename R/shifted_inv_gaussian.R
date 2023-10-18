#' Shifted inverse Gaussian density function
#'
#' @param x Value space of likelihood, x > 0
#' @param mu Mean , mu > 0
#' @param shape Shape, shape unbound
#' @param shift Shift, shift >= 0
#' @param log logical; if TRUE, log(pdf) is returned
#'
#' @return f(x | mu, shape, shift)
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 10, length.out = 1000)
#' plot(x, dshifted_inv_gaussian(x, 1, 1, 1))
dshifted_inv_gaussian <- function(x, mu, shape, shift, log = FALSE) {
  if (isTRUE(any(x <= 0))) {
    stop("Argument has to be > 0")
  }
  if (isTRUE(any(mu <= 0))) {
    stop("Argument has to be > 0")
  }
  if (isTRUE(any(shape <= 0))) {
    stop("Argument has to be > 0")
  }
  if (isTRUE(any(shift < 0))) {
    stop("Shift has to be >= 0")
  }
  brms::dinv_gaussian(x - shift, mu, shape, log)
}

#' Shifted inverse Gaussian RNG function
#'
#' @param n number of observations
#' @param mu Mean, mu > 0
#' @param shape Shape, shape unbound
#' @param shift Shift, shift >= 0
#'
#' @return n samples drawn from the shifted inverse Gaussian distribution
#' @export
#'
#' @examples hist(rshifted_inv_gaussian(100, 1, 1, 1))
rshifted_inv_gaussian <- function(n, mu = 1, shape = 1, shift = 1) {
  if (isTRUE(any(mu <= 0))) {
    stop("Argument has to be > 0")
  }
  if (isTRUE(any(shape <= 0))) {
    stop("Argument has to be > 0")
  }
  if (isTRUE(any(shift < 0))) {
    stop("Shift has to be >= 0")
  }
  brms::rinv_gaussian(n, mu, shape) + shift
}

#' Log-Likelihood of the shifted inverse Gaussian family
#'
#' @param i Indices
#' @param prep brms data
#'
#' @return Log-Likelihood of brms data
log_lik_shifted_inv_gaussian <- function(i, prep) {
  return(
    dshifted_inv_gaussian(
      x = prep$data$Y[i],
      mu = brms::get_dpar(prep, "mu", i = i),
      shape = brms::get_dpar(prep, "shape", i = i),
      shift = brms::get_dpar(prep, "ndt", i = i),
      log = TRUE
    )
  )
}

#' Posterior Prediction of the Shifted Inverse Gaussian family
#'
#' @param i Indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Draws from the posterior predictive distribution
posterior_predict_shifted_inv_gaussian <- function(i, prep, ...) {
  return(
    rshifted_inv_gaussian(
      n = prep$ndraws,
      mu = brms::get_dpar(prep, "mu", i = i),
      shape = brms::get_dpar(prep, "shape", i = i),
      shift = brms::get_dpar(prep, "ndt", i = i)
    )
  )
}

#' Expectations of posterior predictions of the Shifted Inverse Gaussian family
#'
#' @param prep brms data
#'
#' @return Mean of the posterior predictive distribution
posterior_epred_shifted_inv_gaussian <- function(prep) {
  return(brms::get_dpar(prep, "mu") + brms::get_dpar(prep, "ndt"))
}

#' Shifted Inverse Gauss brms family
#'
#' @param link link for mu, default="log"
#' @param link_shape link for the shape, default="log"
#' @param link_ndt link for the shift, default="log"
#'
#' @return Shifted Inverse Gauss brms family
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(
#'   a = a,
#'   y = rshifted_inv_gaussian(
#'     n = 1000, mu = exp(0.5 * a + 1),
#'     shape = 1, shift = 1
#'   )
#' )
#' fit <- brms::brm(
#'   formula = y ~ 1 + a, data = data,
#'   family = shifted_inv_gaussian(), stanvars = shifted_inv_gaussian()$stanvars,
#'   refresh = 0
#' )
#' plot(fit)
shifted_inv_gaussian <- function(link = "log",
                                 link_shape = "log",
                                 link_ndt = "log") {
  family <- brms::custom_family(
    "shifted_inv_gaussian",
    dpars = c("mu", "shape", "ndt"),
    links = c(link, link_shape, link_ndt),
    lb = c(0, 0, 0),
    ub = c(NA, NA, "min_Y"),
    type = "real",
    log_lik = log_lik_shifted_inv_gaussian,
    posterior_predict = posterior_predict_shifted_inv_gaussian,
    posterior_epred = posterior_epred_shifted_inv_gaussian
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real shifted_inv_gaussian_lpdf(real y, real mu, real shape, real ndt) {
          return 0.5 * log(shape / (2 * pi())) - 1.5 * log(y - ndt)
               - 0.5 * shape * square(((y - ndt) - mu) / (mu * sqrt(y - ndt)));
      }",
    block = "functions"
  )
  return(family)
}
