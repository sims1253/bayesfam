# suggested in doi:10.1080/02664763.2018.1511774 for data on the unit interval
# https://www.tandfonline.com/doi/full/10.1080/02664763.2018.1511774


#' Density of Unit Lindley Likelihood
#'
#' @param x Value space of the Unit-Lindley Likelihood family, x e (0, 1)
#' @param mu Mean paraemter of Unit-Lindley Likelihodd, mu e (0, 1)
#' @param log Logical scalar parameter. if TRUE returns log of PDF, default=FALSE
#'
#' @return density f(x | mu)
#' @export
#'
#' @examples x <- seq(from=0.1, to=0.9, length.out=1000)
#' plot(x, dunit_lindley(x, 0.5))
dunit_lindley <- function(x, mu, log=FALSE) {
  if(isTRUE(any(x <= 0 | x >= 1))) {
    stop("The x argument has to be in the interval (0, 1)")
  }
  if(isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mu argument has to be in the interval (0, 1)")
  }

  lpdf <- 2*log1p(-mu)-log(mu)-3*log1p(-x)-x*(1-mu)/(mu*(1-x))
  if(log) {
    return(lpdf)
  }
  else {
    return(exp(lpdf))
  }
}

#' Quantile function of Unit-Lindley Likelihood family
#'
#' @param p Percentile of Unit-Lindley family, p e (0, 1)
#' @param mu Mean parameter, mu e (0, 1)
#'
#' @return QF(p | mu)
#' @export
#'
#' @examples p <- seq(from=0.1, to=0.9, length.out=100)
#' plot(p, qunit_lindley(p, mu=0.5))
qunit_lindley <- function(p, mu) {
  if(isTRUE(any(p <= 0 | p >= 1))) {
    stop("The p argument has to be in the interval (0, 1)")
  }
  if(isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mu argument has to be in the interval (0, 1)")
  }

  lambert_term <- lamW::lambertWm1((1/mu) * (p-1) * exp(-1/mu))

  enumerator <- 1/mu + lambert_term
  divisor <- 1 + lambert_term

  return(enumerator/divisor)
}

#' Unit-Lindley RNG function
#'
#' @param n Number of draws, scalar natural number
#' @param mu Mean argument, mu e (0, 1)
#'
#' @return N number of samples drawn from Unit-Lindley likelihood
#' @export
#'
#' @examples hist(runit_lindley(100, 0.5))
runit_lindley <- function(n, mu) {
  return(qunit_lindley(runif(n), mu))
}

#' Logarithmic Density BRMS vignette of Unitlindley Likelihood
#'
#' @param i Indices of BRMS data
#' @param prep BRMS data
#'
#' @return Log-Likelihood of BRMS data
#' @export
log_lik_unitlindley <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  y <- prep$data$Y[i]
  return(dunit_lindley(y, mu, log=TRUE))
}

#' Posterior mean BRMS vingette of Unitlindley Likelihood
#'
#' @param prep BRMS data
#'
#' @return Mean of the Unitlindley posterior
#' @export
posterior_epred_unitlindley <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

#' Posterior Prediction BRMS Vignette of Unitlindley Distribution
#'
#' @param i Indices of BRMS data
#' @param prep BRMS data
#' @param ... Catchall argument
#'
#' @return  Draws from the Posterior Predictive Distribution
#' @export
posterior_predict_unitlindley <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  return(runit_lindley(prep$ndraws, mu))
}

#' Unitlindley BRMS family
#'
#' @param link link argument of the mean argument, default = log
#'
#' @return Unitlindley BRMS object
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = runit_lindley(n = 1000, mu = inv_logit(0.5 * a + 1)))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'  family = unit_lindley(), stanvars = unit_lindley()$stanvars,
#'  refresh = 0)
#' plot(fit)
unit_lindley <- function(link = "logit") {
  family <- brms::custom_family(
    "unit_lindley",
    dpars = c("mu"),
    links = c(link),
    lb = c(0),
    ub = c(1),
    type = "real",
    log_lik = log_lik_unitlindley,
    posterior_predict = posterior_predict_unitlindley,
    posterior_epred = posterior_epred_unitlindley
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real unit_lindley_lpdf(real y, real mu) {
        return 2*log1m(mu)-log(mu)-3*log1m(y)-y*(1-mu)/(mu*(1-y));
      }
    ",
    block = "functions"
  )
  return(family)
}
