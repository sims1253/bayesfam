#' Probability density function for the Lomax distribution, with Mean parametrization.
#'
#' @param x Model space, defined for x >= 0
#' @param mu Mean parameter of pdf, mu > 0
#' @param alpha Alpha parameter of pdf, alpha > 1
#' @param log Optional log argument, if true, return log(pdf)
#'
#' @details \deqn{f(y | \mu, \alpha) = \frac{\alpha}{\mu(\alpha-1)} (1 + \frac{y}{\mu(\alpha-1)})^{-\alpha-1}}
#'
#' @return PDF of Lomax Distribution
#' @export
#'
#' @examples x <- seq(from = 0, to = 10, length.out = 100)
#' plot(x, dlomax(x, mu = 1, alpha = 2), type = "l")
dlomax <- function(x, mu, alpha, log = FALSE) {
  # check arguments
  if (isTRUE(any(x < 0))) {
    stop("The Lomax PDF is defined only on the positive scale")
  }
  if (isTRUE(alpha <= 1)) {
    stop("The Lomax PDF is defined only for alpha > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("The Lomax PDF is defined only for mu > 0")
  }
  lpdf <- log(alpha) +
    alpha * (log(mu) + log(alpha - 1)) -
    (alpha + 1) * (log(x + (mu * (alpha - 1))))

  # lpdf <- log(alpha) - log(mu) + log(alpha - 1) - (alpha + 1) * log1p(x / (mu * (alpha - 1)))
  if (log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

#' Quantile function for the Lomax distribution, with Mean parametrization.
#'
#' @param p Quantile to be calculated
#' @param mu Median argument of Lomax
#' @param alpha Alpha argument of Gompertz
#'
#' @return Inverse of CDF, calculates a value, given a probability p
#' @export
#'
#' @examples x <- seq(from = 0, to = 1, length.out = 100)
#' plot(x, qlomax(x, mu = 1, alpha = 2), type = "l")
qlomax <- function(p, mu, alpha) {
  # check arguments
  if (isTRUE(any(p < 0))) {
    stop("The Lomax quantile function is defined only on the positive scale")
  }
  if (isTRUE(alpha <= 1)) {
    stop("The Lomax quantile function is defined only for alpha > 1")
  }
  if (isTRUE(mu <= 0)) {
    stop("The Lomax quantile function is defined only for mu > 0")
  }

  # calculate the QDF as inverse of CDF
  return((mu * (alpha - 1)) * ((1 - p)^(-1 / alpha) - 1))
}

#' RNG function for the Lomax distribution, with Mean parametrization.
#'
#' @param n Number of draws
#' @param mu Median argument of Lomax
#' @param alpha Eta argument of Lomax
#'
#' @return A Lomax distributed RNG vector of size n
#' @export
#'
#' @examples hist(log(rlomax(1000, mu = 1, alpha = 2)))
rlomax <- function(n, mu, alpha) {
  # check arguments
  if (isTRUE(mu <= 0)) {
    stop("The Lomax RNG is only defined for mu > 0")
  }
  if (isTRUE(alpha <= 0)) {
    stop("The Lomax RNG is only defined for eta > 0")
  }
  return(qlomax(runif(n, min = 0, max = 1), mu = mu, alpha = alpha))
}

#' Log-Likelihood vignette for the Lomax distribution, with Mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of Lomax given data in prep
log_lik_lomax <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  alpha <- brms::get_dpar(prep, "alpha", i = i)
  y <- prep$data$Y[i]
  return(dlomax(y, mu, alpha, log = TRUE))
}

#' Posterior-Prediction vignette for the Lomax distribution, with Mean parametrization.
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ...
#'
#' @return Posterior prediction of Lomax, given data in prep
posterior_predict_lomax <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  alpha <- brms::get_dpar(prep, "alpha", i = i)
  return(rgompertz(prep$ndraws, mu, alpha))
}

#' Expectation-Predict vignette for the Lomax distribution, with Mean parametrization.
#'
#' @param prep BRMS data
#'
#' @return Recover the given mean of data prep
posterior_epred_lomax <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  return(mu)
}

#' Lomax Stan-implementation in Mean parametrization.
#'
#' @param link Link function for function
#' @param link_alpha Link function for eta argument
#'
#' @return BRMS Lomax distribution family
#' @export
#'
#' @examples # Running the example might take a while and may make RStudio unresponsive.
#' # Just relax and grab a cup of coffe or tea in the meantime.
#' a <- rnorm(1000)
#' data <- list(a = a, y = rlomax(1000, exp(0.5 * a + 1), 2))
#' # refresh = 0 supresses chain updates
#' fit1 <- brms::brm(y ~ 1 + a, data = data,
#'  family = lomax(), stanvars = lomax()$stanvars,
#'  refresh = 0)
#' plot(fit1)
lomax <- function(link = "log", link_alpha = "log1p") {
  family <- brms::custom_family(
    "lomax",
    dpars = c("mu", "alpha"),
    links = c(link, link_alpha),
    lb = c(0, 1),
    ub = c(NA, NA),
    type = "real",
    log_lik = log_lik_lomax,
    posterior_predict = posterior_predict_lomax,
    posterior_epred = posterior_epred_lomax
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real lomax_lpdf(real y, real mu, real alpha) {
        return(log(alpha) +
               alpha * (log(mu) + log(alpha - 1)) -
               (alpha + 1) * (log(y + (mu * (alpha - 1)))));
      }

      real lomax_rng(real mu, real alpha) {
        return ((mu * (alpha - 1)) * ((1 - uniform_rng(0, 1))^(-1/alpha) - 1));
      }",
    block = "functions"
  )
  return(family)
}
