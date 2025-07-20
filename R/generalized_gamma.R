#' Generalized gamma distribution
#'
#' @source Bases on flexsurv
#' (<https://github.com/chjackson/flexsurv/tree/master>) by Christopher
#' Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>. Inspired by a blog post by
#' Demetri Pananos (\url{https://dpananos.github.io/posts/2023-12-02-gen-gamma/})
#' and code by Krzysztof Sakrejda (<https://github.com/sakrejda/tooling>).
#'
#' @param x Value, x > 0.
#' @param n number of observations.
#' @param mu Vector of ``location'' parameters.
#' @param sigma Vector of ``scale'' parameters.  Note the inconsistent
#' meanings of the term ``scale'' - this parameter is analogous to the
#' (log-scale) standard deviation of the log-normal distribution, ``sdlog'' in
#' [dlnorm()], rather than the ``scale'' parameter of the gamma
#' distribution [dgamma()]. Constrained to be positive.
#' @param Q Vector of shape parameters.
#' @param log logical; if TRUE the log-pdf is returned
#' @param link Link function for mu
#' @param link_sigma Link function for sigma
#' @param link_Q Link function for Q
#'
#' @return `dgeneralized_gamma` gives the density,
#' `rgeneralized_gamma` generates random deviates,
#'
#' @references Prentice, R. L. (1974). A log gamma model and its maximum
#' likelihood estimation. Biometrika 61(3):539-544.
#'
#' Farewell, V. T. and Prentice, R. L. (1977). A study of
#' distributional shape in life testing. Technometrics 19(1):69-75.
#'
#' Lawless, J. F. (1980). Inference in the generalized gamma and log
#' gamma distributions.  Technometrics 22(3):409-419.
#'
#' Cox, C., Chu, H., Schneider, M. F. and Muñoz, A. (2007).
#' Parametric survival analysis and taxonomy of hazard functions for
#' the generalized gamma distribution.  Statistics in Medicine
#' 26:4252-4374
#'
#' Stacy, E. W. (1962). A generalization of the gamma distribution.
#' Annals of Mathematical Statistics 33:1187-92
#' @name generalized_gamma
NULL

##' @export
##' @rdname generalized_gamma
dgeneralized_gamma <- function(x, mu = 0, sigma = 1, Q, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("The generalized gamma distribution is only defined for x > 0")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("The generalized gamma distribution is only defined for sigma > 0")
  }

  # w = (log(x) - mu) / sigma
  # lpdf = (k - 0.5) * log(k) - log(sigma) - lgamma(k) +
  #        (sqrt(k) * w - k * exp(1 / sqrt(k) * w)) - log(x)

  if (Q != 0) {
    qi <- 1 / (Q * Q)
    qw <- Q * ((log(x) - mu) / sigma)
    lpdf <- -log(sigma * x) +
      log(abs(Q)) * (1 - 2 * qi) +
      qi * (qw - exp(qw)) -
      lgamma(qi)
  } else {
    lpdf <- dlnorm(x, mu, sigma, 1)
  }

  # return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' @export
#' @rdname generalized_gamma
rgeneralized_gamma <- function(n, mu = 0, sigma = 1, Q) {
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sigma) == 1) {
    sigma <- rep(sigma, n)
  }
  if (length(Q) == 1) {
    Q <- rep(Q, n)
  }
  if (any(sigma <= 0)) {
    stop(
      "The generalized gamma distribution is only defined for positive values
         of sigma!"
    )
  }

  q0 <- Q == 0
  out <- vector(mode = "numeric", length = n)
  out[which(q0)] <- rlnorm(sum(q0), mu[which(q0)], sigma[which(q0)])

  out[which(!q0)] <- exp(
    mu[which(!q0)] +
      sigma[which(!q0)] *
        (log(
          Q[which(!q0)]^2 *
            rgamma(
              sum(!q0),
              1 / Q[which(!q0)]^2,
              1
            )
        ) /
          Q[which(!q0)])
  )

  # This sometimes results in zeroes, which we don't want, so we truncate at
  # machine precision

  out[which(out == 0)] <- .Machine$double.xmin * 1e8
  return(out)
}

#' Log-Likelihood of the generalized gamma distribution
#'
#' @param i brms indices
#' @param prep brms data
#'
#' @return Log-Likelihood of generalized gamma given data in prep
log_lik_generalized_gamma <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  Q <- brms::get_dpar(prep, "Q", i = i)
  y <- prep$data$Y[i]
  return(dgeneralized_gamma(x = y, mu = mu, sigma = sigma, Q = Q, log = TRUE))
}

#' posterior_predict for the generalized gamma distribution
#'
#' @param i brms indices
#' @param prep brms data
#' @param ... Catchall argument
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_generalized_gamma <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  Q <- brms::get_dpar(prep, "Q", i = i)
  return(rgeneralized_gamma(n = prep$ndraws, mu = mu, sigma = sigma, Q = Q))
}

#' posterior_epred for the generalized gamma distribution
#'
#' @param prep brms data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_generalized_gamma <- function(prep) {
  stop("There is no implementation of posterior_epred!")
}


#' @export
#' @rdname generalized_gamma
generalized_gamma <- function(
  link = "log",
  link_sigma = "log",
  link_Q = "log"
) {
  family <- brms::custom_family(
    "generalized_gamma",
    dpars = c("mu", "sigma", "Q"),
    links = c(link, link_sigma, link_Q),
    lb = c(-NA, 0, -NA),
    ub = c(NA, NA, NA),
    type = "real",
    log_lik = log_lik_generalized_gamma,
    posterior_predict = posterior_predict_generalized_gamma,
    posterior_epred = posterior_epred_generalized_gamma
  )
  family$stanvars <- brms::stanvar(
    scode = "
    real generalized_gamma_lpdf(real y, real mu, real sigma, real Q) {
      if (Q!=0) {
        real qi = 1/(Q * Q);
        real qw = Q * ((log(y) - mu) / sigma);
        return -log(sigma*y) +
          log(abs(Q)) * (1 - 2 * qi) +
          qi * (qw - exp(qw)) - lgamma(qi);
      } else {
        return lognormal_lpdf(y| mu, sigma);
      }
  }",
    block = "functions"
  )
  return(family)
}
