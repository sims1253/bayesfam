# Euler-Mascheroni constant used in a few equations in Gumbel
euler_mascheroni <- 0.57721566490153

#' Probability density function for the gumbel distribution
#'
#' @details The gumbel distribution has density
#' \deqn{f(y | \mu, \sigma) = \frac{1}{\sigma} exp(-(z + e^{-z}))}
#' @details Where z is the linear transformation
#' \deqn{z(y, \mu, \sigma) = \frac{y - \mu}{\sigma} + \gamma}
#' @details Where Gamma refers to the Euler-Mascheroni constant
#'
#' @param x Value, unbound
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#' @param log Optional argument. If true, returns the log density.
#'
#' @return Density of the pdf given x, mu and sigma
#' @export
#'
#' @examples x <- seq(from = -5, to = 10, length.out = 1000)
#' plot(x, dgumbel_mean(x, mu = 2, sigma = 2), type = "l")
dgumbel_mean <- function(x, mu, sigma, log = FALSE) {
  # check the arguments
  if (isTRUE(any(sigma <= 0))) {
    stop("gumbel is only defined for sigma > 0")
  }
  # Maybe overkill?
  if (!lenEqual(list_of_vectors = list(x, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("Gumbel argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }
  if (!isLogic_len(log)) {
    stop("the log argument of a density has to be a scalar boolean")
  }

  # location <- mu - sigma * euler_mascheroni
  z <- (x - mu) / sigma + euler_mascheroni
  lpdf <- -log(sigma) - (z + exp(-z))

  #return either the log or the pdf itself, given the log-value
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function of the gumbel distribution
#'
#' @param p quantile value, 0 < p < 1
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#'
#' @return Quantiles of the gumbel distribution
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 100)
#' plot(x, qgumbel_mean(x, mu = 2, sigma = 2), type = "l")
qgumbel_mean <- function(p, mu, sigma) {
  # check the arguments
  if (isTRUE(any(p <= 0 | p >= 1))) {
    stop("quantile value has to be between 0 and 1")
  }
  if (isTRUE(any(sigma <= 0))) {
    stop("gumbel is only defined for sigma > 0")
  }
  if (!lenEqual(list_of_vectors = list(p, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
    stop("Gumbel argument vectors could not be matched. May be due to wrong type,
         or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
  }

  location <- mu - sigma * euler_mascheroni
  return(location - sigma * log(-log(p)))
}

#' RNG for the gumbel distribution
#'
#' @param n Number of samples
#' @param mu Mean, unbound
#' @param sigma Scale, sigma > 0
#'
#' @return Random numbers from the gumbel distribution.
#' @export
#'
#' @examples hist(rgumbel_mean(100, mu = 2, sigma = 2))
rgumbel_mean <- function(n, mu = 0, sigma = 1) {
  # check the arguments
  if (!isNat_len(n)) {
    stop("The number RNG-samples has to be a scalar natural")
  }
  return(qgumbel_mean(runif(n, min = 0, max = 1), mu, sigma))
}

#' Log-Likelihood of the gumbel distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of gumbel given data in prep
log_lik_gumbel_mean <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dgumbel_mean(y, mu, sigma, log = TRUE))
}


#' posterior_predict for the gumbel distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ... catchall argument
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_gumbel_mean <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rgumbel_mean(prep$ndraws, mu, sigma))
}

#' posterior_epred for the gumbel distribution
#'
#' @param prep BRMS data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_gumbel_mean <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' Gumbel brms custom family, custom implementation vs Stan
#'
#' @param link Link function for function
#' @param link_sigma Link function for sigma argument
#'
#' @return BRMS Gumbel distribution family
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = rgumbel_mean(n = 1000, mu = a + 2, sigma = 2))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'   family = gumbel_mean(), stanvars = gumbel_mean()$stanvars,
#'   refresh = 0)
#' plot(fit)
gumbel_mean <- function(link = "identity", link_sigma = "log") {
  family <- brms::custom_family(
    "gumbel_mean",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(-NA, 0),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_gumbel_mean,
    posterior_predict = posterior_predict_gumbel_mean,
    posterior_epred = posterior_epred_gumbel_mean
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real gumbel_mean_lpdf(real y, real mu, real sigma) {
        real loc = mu - sigma * 0.57721566490153; // Euler-Mascheroni
        return gumbel_lpdf(y | loc, sigma);
      }

      real gumbel_mean_rng(real mu, real sigma) {
        real loc = mu - sigma * 0.57721566490153;
        return gumbel_rng(loc, sigma);
      }",
    block = "functions"
  )
  return(family)
}
