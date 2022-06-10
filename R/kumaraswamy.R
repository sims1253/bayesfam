#' Title
#'
#' @param x
#' @param mu
#' @param p
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dkumaraswamy <- function(x, mu, p, log = FALSE) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop("The mean must be in (0,1).")
  }
  mu[which(mu > 0.999999)] <- 0.999999
  mu[which(mu < 0.000001)] <- 0.000001
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  logpdf <- log(p) +
    log(log(2)) -
    log(-(log1p(-mu^p))) +
    (p - 1) * log(x) +
    ((-(log(2) / log1p(-mu^p))) - 1) * log1p(-x^p)

  if (log) {
    return(logpdf)
  } else {
    return(exp(logpdf))
  }
}

#' Title
#'
#' @param n
#' @param mu
#' @param p
#'
#' @return
#' @export
#'
#' @examples
rkumaraswamy <- function(n, mu, p) {
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop("The mean must be in (0,1).")
  }
  mu[which(mu > 0.999999)] <- 0.999999
  mu[which(mu < 0.000001)] <- 0.000001
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  q <- -(log(2) / log1p(-mu^p))
  return(
    (1 - (1 - runif(n, min = 0, max = 1))^(1 / q))^(1 / p)
  )
}

#' Title
#'
#' @param u
#' @param mu
#' @param p
#'
#' @return
#' @export
#'
#' @examples
qkumaraswamy <- function(u, mu = 0.5, p = 1) {
  if (isTRUE(any(u <= 0 | u >= 1))) {
    stop("u must be in (0,1).")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  return(
    (1 -
      (1 - u)^(1 /
        (-(log(2) / log1p(-mu^p))))
    )^(1 / p)
  )
}

#' Title
#'
#' @param x
#' @param mu
#' @param p
#'
#' @return
#' @export
#'
#' @examples
pkumaraswamy <- function(x, mu = 0.5, p = 1) {
  if (isTRUE(any(x <= 0 | x >= 1))) {
    stop("x must be in (0,1).")
  }
  if (isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The median must be in (0,1).")
  }
  if (isTRUE(any(p <= 0))) {
    stop("P must be above 0.")
  }
  q <- -(log(2) / log1p(-mu^p))
  return(1 + (x^p - 1)^q)
}

#' Title
#'
#' @param i
#' @param prep
#'
#' @return
#'
#' @examples
log_lik_kumaraswamy <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  y <- prep$data$Y[i]
  return(dkumaraswamy(y, mu, p, log = TRUE))
}

#' Title
#'
#' @param i
#' @param prep
#' @param ...
#'
#' @return
#'
#' @examples
posterior_predict_kumaraswamy <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  p <- brms::get_dpar(prep, "p", i = i)
  return(rkumaraswamy(prep$ndraws, mu, p))
}

#' Title
#'
#' @param prep
#'
#' @return
#'
#' @examples
posterior_epred_kumaraswamy <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  p <- brms::get_dpar(prep, "p")
  q <- -(log(2) / log1p(-mu^p))
  return(q * beta((1 + 1 / p), q))
}

#' Title
#'
#' @param link
#' @param link_p
#'
#' @return
#' @export
#'
#' @examples
kumaraswamy <- function(link = "logit", link_p = "log") {
  family <- brms::custom_family(
    "kumaraswamy",
    dpars = c("mu", "p"),
    links = c(link, link_p),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "real",
    log_lik = log_lik_kumaraswamy,
    posterior_predict = posterior_predict_kumaraswamy,
    posterior_epred = posterior_epred_kumaraswamy
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real kumaraswamy_lpdf(real y, real mu, real p) {
         return  (log(p) + log(log(2)) - log(-(log1m(mu^p))) + (p-1) * log(y) +
                 ((-(log(2)/log1m(mu^p)))-1) * log1m(y^p));
      }

      real kumaraswamy_rng(real mu, real p) {
         return ((1-(1-uniform_rng(0, 1))^(1/(-(log(2)/log1m(mu^p)))))^(1/p));
      }",
    block = "functions"
  )
  return(family)
}
