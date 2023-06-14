
#' Density of Inverse Gauss Likelihood with shift parameter.
#'
#' @param x Value space of likelihood, x > 0
#' @param mu Mean parameter? of likelihod function, mu > 0
#' @param shape Shape parameter of likelihood, shape unbound
#' @param shift Shift paramter of likelihood, shift <= 0
#' @param log Logical paramter, if TRUE returns log PDF, default=FALSE
#'
#' @return f(x | mu, shape, shift)
#' @export
#'
#' @examples x <- seq(from=0.1, to=10, length.out=1000)
#' plot(x, dshifted_inv_gaussian(x, 1, 1, 1))
dshifted_inv_gaussian <- function(x, mu, shape, shift, log = FALSE) {
  if(isTRUE(any(x <= 0))) {
    stop("Argument has to be > 0")
  }
  if(isTRUE(any(mu <= 0))) {
    stop("Argument has to be > 0")
  }
  if(isTRUE(any(shape <= 0))) {
    stop("Argument has to be > 0")
  }
  # if(isTRUE(any(shift > 0))) {
  #   stop("Shift has to be <= 0")
  # }
  brms::dinv_gaussian(x-shift, mu, shape, log)
}

#' RNG function of the Shifted Inverse Gauss Likelihood
#'
#' @param n Number samples to draw, natural scalar
#' @param mu Mean parameter? of likelihod function, mu > 0
#' @param shape Shape parameter of likelihood, shape unbound
#' @param shift Shift paramter of likelihood, shift <= 0
#'
#' @return N samples drawn of shifted inverse Gauss
#' @export
#'
#' @examples hist(rshifted_inv_gaussian(100, 1, 1, -1))
rshifted_inv_gaussian <- function(n, mu = 1, shape = 1, shift = 0) {
  if(isTRUE(any(mu <= 0))) {
    stop("Argument has to be > 0")
  }
  if(isTRUE(any(shape <= 0))) {
    stop("Argument has to be > 0")
  }
  # if(isTRUE(any(shift > 0))) {
  #   stop("Shift has to be <= 0")
  # }
  brms::rinv_gaussian(n, mu, shape) + shift
}

#' Posterior mean BRMS vingette of Shifted Inverse Gauss Likelihood
#'
#' @param prep BRMS data
#'
#' @return Mean of the shifted inverse Gauss posterior
#' @export
posterior_epred_shifted_inv_gaussian <- function(prep) {
    #with(prep$dpars, mu - ndt)
  mu <- brms::get_dpar(prep, "mu", i = i)
  shift <- brms::get_dpar(prep, "ndt", i = i)
  return(mu + shift)
}

#' Posterior Prediction BRMS Vignette of Shifted Inverse Gauss
#'
#' @param i Indices of BRMS data
#' @param prep BRMS data
#' @param ... Catchall argument
#'
#' @return  Draws from the Posterior Predictive Distribution
#' @export
posterior_predict_shifted_inv_gaussian <- function(i, prep, ...) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    return(rshifted_inv_gaussian(prep$ndraws, mu, shape, ndt))
}

#' Logarithmic Density BRMS vignette of shifted inverse Gauss Likelihood
#'
#' @param i Indices of BRMS data
#' @param prep BRMS data
#'
#' @return Log-Likelihood of BRMS data
#' @export
log_lik_shifted_inv_gaussian <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    y <- prep$data$Y[i]
    return(dshifted_inv_gaussian(y, mu, shape, ndt, log = TRUE))
}

#' Shifted Inverse Gauss BRMS function
#'
#' @param link Central paramter link function, default="log"
#' @param link_shape Shape parameter link function, default="log"
#' @param link_ndt Shift parameter link function, default="log"
#'
#' @return Shifted Inverse Gauss BRMS family
#' @export
#'
#' @examples a <- rnorm(1000)
#' data <- list(a = a, y = rshifted_inv_gaussian(n=1000, mu=exp(0.5*a + 1), shape=1, shift=1))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'  family = shifted_inv_gaussian(), stanvars = shifted_inv_gaussian()$stanvars,
#'  refresh = 0)
#' plot(fit)
shifted_inv_gaussian <- function(link = "log", link_shape = "log", link_ndt = "log"){
  # test as in BRMS additional families
  family <- brms::custom_family(
      "shifted_inv_gaussian",
      dpars = c("mu", "shape", "ndt"),
      links = c(link, link_shape, link_ndt),
      lb = c(0, 0, 0),
      ub = c(NA, NA, NA),
      type = "real",
      log_lik = log_lik_shifted_inv_gaussian,
      posterior_predict = posterior_predict_shifted_inv_gaussian,
      posterior_epred = posterior_epred_shifted_inv_gaussian
    )
  family$stanvars <- brms::stanvar(
    scode = "
      real shifted_inv_gaussian_lpdf(real y, real mu, real shape, real ndt) {

          real x = y - ndt;
          return 0.5 * log(shape / (2 * pi())) - 1.5 * log(x)
               - 0.5 * shape * square((x - mu) / (mu * sqrt(x)));
      }

      real shifted_inv_gaussian_rng(real mu, real shape, real ndt) {
          real y = normal_rng(0, 1) ^ 2;
          real x = mu + (mu^2 * y) / (2 * shape) - mu / (2 * shape) *
                    sqrt(4 * mu * shape * y + mu^2 * y^2);
          real z = uniform_rng(0,1);

          if(z <= mu / (mu + x))
            return x + ndt;
          else
            return mu^2 / x + ndt;
      }",
    block = "functions"
  )
  return(family)
}
