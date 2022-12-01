#' Simplex density function in mean parametrisation.
#'
#' @param x value space, x e (0, 1)
#' @param mu Median parameter of pdf, mu e (0, 1)
#' @param sigma shape parameter, sigma unbound
#' @param log if true, returns log(pdf). Normally FALSE.
#'
#' @details \deqn{f(y) = (2 \pi \sigma^2(y(1-y))^3)^{-\frac{1}{2}} exp(-(\frac{y-\mu}{\mu(1-\mu)})^2 \frac{1}{2y(1-y)\sigma^2} )}
#'
#' @return f(x | mu, sigma)
#' @export
#'
#' @examples x <- seq(from = 0.01, to = 0.99, length.out = 1000)
#' plot(x, dsimplex(x, mu = 0.7, sigma = 2), type = "l")
dsimplex <- function(x, mu, sigma, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop("The mean must be in (0,1).")
  }
  mu[which(mu > 0.999999)] <- 0.999999
  mu[which(mu < 0.000001)] <- 0.000001
  if (isTRUE(any(sigma == 0))) {
    stop("sigma can not be 0.")
  }
  result <- (-0.5) * (
    log(2) +
      log(pi) +
      2 * log(sigma) +
      3 * (log(x) + log1p(-x))
  ) +
    ((-1 / (2 * sigma^2)) * (
      ((x - mu)^2) /
        (x * (1 - x) * mu^2 * (1 - mu)^2)
    ))
  if (log) {
    return(result)
  } else {
    return(exp(result))
  }
}


#' RNG for the inverse gaussian distribution
#'
#' Based on code from simplexreg
#' Peng Zhang, Zhenguo Qiu, Chengchun Shi (2016). simplexreg: An R
#' Package for Regression Analysis of Proportional Data Using the
#' Simplex Distribution. Journal of Statistical Software, 71(11), 1-21.
#' doi:10.18637/jss.v071.i11
#'
#' @param epsilon epsilon parameter of the inverse gaussian distribution
#' @param Tau tau parameter of the inverse gaussian distribution
#'
#' @return a single sample from the specified inverse gaussian
rIG <-
  function(epsilon, Tau) {
    ## generating random number from inverse-gaussian dist'n
    z <- rchisq(1, 1)
    ss <- sqrt(4 * epsilon * z / Tau + (epsilon * z)^2)
    z1 <- epsilon + (epsilon^2) * Tau * z / 2 - (epsilon * Tau / 2) * ss
    u1 <- runif(1, 0, 1)
    xxx <- z1
    if (u1 > (epsilon / (epsilon + z1))) {
      xxx <- (epsilon^2) / z1
    }
    return(as.numeric(xxx))
  }


#' RNG for a mixture of inverse gaussian distributions
#'
#' Based on code from simplexreg
#' Peng Zhang, Zhenguo Qiu, Chengchun Shi (2016). simplexreg: An R
#' Package for Regression Analysis of Proportional Data Using the
#' Simplex Distribution. Journal of Statistical Software, 71(11), 1-21.
#' doi:10.18637/jss.v071.i11
#'
#' @param epsilon epsilon parameters of the mixture parts
#' @param Tau tau parameters of the mixture parts
#' @param mu mu parameters of the mixture parts
#'
#' @return a single sample for the specified mixture
rMIG <-
  function(epsilon, Tau, mu) {
    ## generating random number from inverse-gaussian mixture dist'n
    x1 <- rIG(epsilon, Tau)
    x2 <- rchisq(1, 1)
    x3 <- x2 * Tau * (epsilon^2)
    u2 <- runif(1, 0, 1)
    xx <- x1
    if (u2 < mu) {
      xx <- x1 + x3
    }
    return(as.numeric(xx))
  }


#' Simplex RNG function in Median parametrization.
#'
#' Based on code from simplexreg
#' Peng Zhang, Zhenguo Qiu, Chengchun Shi (2016). simplexreg: An R
#' Package for Regression Analysis of Proportional Data Using the
#' Simplex Distribution. Journal of Statistical Software, 71(11), 1-21.
#' doi:10.18637/jss.v071.i11
#'
#' @param n Number of samples to draw, as a natural number scalar.
#' @param mu Mean parameter, mu e (0, 1)
#' @param sigma shape parameter, Sigma unbound
#'
#' @return n samples in Simplex distribution.
#' @export
#'
#' @examples hist(rsimplex(10000, mu = 0.7, sigma = 2))
rsimplex <-
  function(n, mu, sigma) {
    ## generating random number from simplex dist'n
    ## by transformation from inverse-gaussian mixture dist'n
    if (any(mu <= 0 | mu >= 1)) {
      stop("The mean must be in (0,1).")
    }
    mu[which(mu > 0.999999)] <- 0.999999
    mu[which(mu < 0.000001)] <- 0.000001
    if (any(sigma == 0)) {
      stop("sigma can not be 0.")
    }

    if (length(mu) == 1) {
      mu <- rep(mu, n)
    }
    if (length(sigma) == 1) {
      sigma <- rep(sigma, n)
    }
    epsilon <- mu / (1 - mu)
    Tau <- sigma^2 * ((1 - mu)^2)
    yy <- rep(0, n)
    for (i in 1:n) {
      x <- rMIG(epsilon[i], Tau[i], mu[i])
      yy[i] <- x / (1 + x)
    }
    return(as.vector(yy))
  }


#' Posterior prediction vignette for the Simplex distribution, in Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of Simplex, given data in prep
posterior_predict_simplex <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rsimplex(prep$ndraws, mu, sigma))
}


#' Log-Likelihood vignette for the Simplex distribution, in Median parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of Simplex given data in prep
log_lik_simplex <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dsimplex(y, mu, sigma, log = TRUE))
}


#' Posterior expected value prediction of the Simplex implementation.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
posterior_epred_simplex <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' Simplex BRMS-implementation in median parametrization.
#'
#' @param link Link function for function
#' @param link_sigma Link function for sigma argument
#'
#' @return BRMS Beta-Custom distribution family
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' a <- rnorm(1000)
#' data <- list(a = a, y = rsimplex(1000, brms::inv_logit_scaled(0.5 * a + 1), 2))
#' # refresh = 0 supresses chain updates
#' fit1 <- brms::brm(y ~ 1 + a, data = data,
#'  family = simplex(), stanvars = simplex()$stanvars,
#'  refresh = 0)
#' plot(fit1)
simplex <- function(link = "logit", link_sigma = "identity") {
  family <- brms::custom_family(
    "simplex",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, -NA),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_simplex,
    posterior_predict = posterior_predict_simplex,
    posterior_epred = posterior_epred_simplex
  )
  family$stanvars <- brms::stanvar(
    scode = "
        real simplex_lpdf(real y, real mu, real sigma) {
           return (-0.5) * (
                     log(2) +
                     log(pi()) +
                     2 * log(sigma) +
                     3 * (log(y) + log1m(y))
                  ) + (
                     (-1 / (2 * sigma^2)) *
                     (
                        ((y-mu)^2) /
                        (y * (1-y) * mu^2 * (1-mu)^2)
                     )
                  );
        }",
    block = "functions"
  )
  return(family)
}
