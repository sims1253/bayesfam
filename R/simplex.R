#' Title
#'
#' @param x
#' @param mu
#' @param sigma
#' @param log
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param epsilon
#' @param Tau
#'
#' @return
#'
#' @examples
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


#' Title
#'
#' @param epsilon
#' @param Tau
#' @param mu
#'
#' @return
#'
#' @examples
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


#' Title
#'
#' @param n
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
rsimplex <-
  function(n, mu, sigma) {
    ## generating random number from simplex dist'n
    ## by transformation from inverse-gaussian mixture dist'n
    if (any(mu < 0 | mu > 1)) {
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


#' Title
#'
#' @param i
#' @param prep
#' @param ...
#'
#' @return
#'
#' @examples
posterior_predict_simplex <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  return(rsimplex(prep$ndraws, mu, sigma))
}


#' Title
#'
#' @param i
#' @param prep
#'
#' @return
#'
#' @examples
log_lik_simplex <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  y <- prep$data$Y[i]
  return(dsimplex(y, mu, sigma, log = TRUE))
}


#' Title
#'
#' @param prep
#'
#' @return
#'
#' @examples
posterior_epred_simplex <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' Title
#'
#' @param link
#' @param link_sigma
#'
#' @return
#' @export
#'
#' @examples
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
