#' Logit link function
#'
#' @param x value of x to be transformed, x e (0, 1)
#'
#' @return logit value of x, x unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, logit(x), type = "l")
logit <- function(x) {
  if (any(x < 0 | x > 1)) {
    stop("The logit link is only defined between 0 and 1!")
  }
  return(qlogis(x))
}

#' Logit response function
#'
#' @param x value of x to be transformed, x unbound
#'
#' @return inverse logit value of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -5, to = 5, length.out = 100)
#' plot(x, inv_logit(x), type = "l")
inv_logit <- function(x) {
  return(plogis(x))
}

#' Complementary-Log-Log link function
#'
#' @param x value of x to be transformed, x e (0, 1)
#'
#' @return cloglog value of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, cloglog(x), type = "l")
cloglog <- function(x) {
  if (any(x < 0 | x > 1)) {
    stop("The logit link is only defined between 0 and 1!")
  }
  log(-log1p(-x))
}


#' Complementary-Log-Log response function
#'
#' @param x value of x to be transformed, x unbound
#'
#' @return inverse-cloglog value of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -3, to = 1, length.out = 100)
#' plot(x, inv_cloglog(x), type = "l")
inv_cloglog <- function(x) {
  return(1 - exp(-exp(x)))
}


#' Cauchit link function
#'
#' @param x value of x to be transformed, not defined for x of whole integer
#'
#' @return cauchit value of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 0.9, length.out = 100)
#' plot(x, cauchit(x), type = "l")
cauchit <- function(x) {
  if (any(x < 0 | x > 1)) {
    stop("The logit link is only defined between 0 and 1!")
  }
  return(qcauchy(x))
}


#' Cauchit response function
#'
#' @param x value of x to be transformed, any real scalar or vector allowed
#'
#' @return inverse cauchit of x, result is e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -10, to = 10, length.out = 100)
#' plot(x, inv_cauchit(x), type = "l")
inv_cauchit <- function(x) {
  return(pcauchy(x))
}

#' Gaussion Error function
#'
#' @param x value to be transformed, x unbound
#'
#' @return erf function of x, result e (0, 1)
#' @export
#'
#' @examples x <- seq(from = -2, to = 2, length.out = 100)
#' plot(x, erf(x), type = "l")
erf <- function(x) {
  return(2 * pnorm(x * sqrt(2)) - 1)
}

inv_erf <- function(x) {
  return(qnorm((x + 1) / 2) / sqrt(2))
}

#' Softplus link function
#'
#' @param x value to be transformed, x positive unbound
#'
#' @return softplus function of x, result is unbound
#' @export
#'
#' @examples x <- seq(from = 0, to = 5, length.out = 100)
#' plot(x, softplus(x), type = "l")
softplus <- function(x) {
  return(x + log(-expm1(-x)))
}


#' Softplus response function
#'
#' @param x value to be transformed, x is unbound
#'
#' @return inv_softplus of x, result is positive unbound
#' @export
#'
#' @examples x <- seq(from = 0.1, to = 5, length.out = 100)
#' plot(x, softplus(x), type = "l")
inv_softplus <- function(x) {
  return(pmax(x, 0) + log1p(exp(-abs(x))))
}

#' Symlog link function
#'
#' @source Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
#'         Mastering Diverse Domains through World Models.
#'         (<https://doi.org/10.48550/arXiv.2301.04104>)
#'
#' @param x value to be transformed, x is unbound
#'
#' @return symlog of x, result is unbound
#' @export
#'
#' @examples
#' symlog(0)
#' symlog(1e10)
#' symlog(-1e10)
symlog <- function(x) {
  return(sign(x) * log1p(abs(x)))
}


#' Symlog response function
#'
#' @source Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
#'         Mastering Diverse Domains through World Models.
#'         (<https://doi.org/10.48550/arXiv.2301.04104>)
#'
#' @param x value to be transformed, x is unbound
#'
#' @return inv_symlog of x, result is unbound
#' @export
#'
#' @examples
#' inv_symlog(0)
#' inv_symlog(10)
#' inv_symlog(-10)
inv_symlog <- function(x) {
  return(sign(x) * (exp(abs(x)) - 1))
}

#' Logarithm of the sum of exponentials.
#'
#' A more numerically stable equivalent to `log(sum(exp(x)))`
#'
#' @details An empty vector returns `-Inf` (the sum is empty). A non-finite
#'   maximum is returned as is: `NA`/`NaN` entries propagate, otherwise any
#'   `Inf` entry gives `Inf` and all-`-Inf` input gives `-Inf`.
#'
#' @source https://en.wikipedia.org/wiki/LogSumExp#log-sum-exp_trick_for_log-domain_calculations
#' @param x a vector of values
#' @return log(sum(exp(x)))
#' @export
logsumexp <- function(x) {
  if (!length(x)) {
    return(-Inf)
  }
  y <- max(x)
  if (!is.finite(y)) {
    return(y)
  }
  return(y + log(sum(exp(x - y))))
}

#' Mean of a transformed normal variable.
#'
#' Computes E[inv_link(Z)] for Z ~ N(mu, sigma) with trapezoid sums on fixed
#' nodes, vectorized over `mu` and `sigma`. This is the expectation that
#' `posterior_epred` promises for the families that transform a normal
#' variable (cauchitnormal, cloglognormal, softplusnormal).
#'
#' For sigma <= 1 it integrates inv_link(mu + sigma * t) over t ~ N(0, 1) with
#' step 0.25. A larger sigma makes that integrand steep, so then it uses
#' `dual` if given and otherwise halves the step for every doubling of sigma.
#' `dual` describes the same expectation as a sum over nodes `l` with
#' `weights` of `kernel(mu, sigma, l)`: if inv_link is the CDF F of some L,
#' E[F(Z)] = P(L < Z) = E[pnorm((mu - L) / sigma)], which is smooth in L.
#'
#' @param mu Location of the normal variable, unbound
#' @param sigma Scale of the normal variable, sigma >= 0
#' @param inv_link Vectorized response function
#' @param dual NULL or a list of `nodes`, `weights` and `kernel`
#'
#' @return E[inv_link(Z)], same dimensions as `mu`.
#' @noRd
latent_normal_mean <- function(mu, sigma, inv_link, dual = NULL) {
  sigma <- sigma + 0 * mu
  mu <- mu + 0 * sigma
  out <- 0 * mu
  wide <- !is.null(dual) & sigma > 1
  if (any(wide)) {
    m <- mu[wide]
    s <- sigma[wide]
    total <- 0
    for (j in seq_along(dual$nodes)) {
      total <- total + dual$weights[j] * dual$kernel(m, s, dual$nodes[j])
    }
    out[wide] <- total
  }
  # sum the narrow draws separately for each step size
  level <- pmax(0, ceiling(log2(sigma)))
  for (k in unique(level[!wide])) {
    narrow <- !wide & level == k
    h <- 0.25 / 2^k
    m <- mu[narrow]
    s <- sigma[narrow]
    total <- 0
    for (t in seq(-9, 9, by = h)) {
      total <- total + h * dnorm(t) * inv_link(m + s * t)
    }
    out[narrow] <- total
  }
  return(out)
}
