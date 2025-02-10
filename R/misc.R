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
  return(qnorm(x / 2 + 1) / sqrt(2))
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
  return(log(exp(x) - 1))
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
  return(log(exp(x) + 1))
}

#' Symlog link function
#'
#' @source Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
#'         Mastering Diverse Domains through World Models.
#'         (\url{https://doi.org/10.48550/arXiv.2301.04104})
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
#'         (\url{https://doi.org/10.48550/arXiv.2301.04104})
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
#' A more numerically stable equivalent to \code{log(sum(exp(x)))}
#'
#' @source https://en.wikipedia.org/wiki/LogSumExp#log-sum-exp_trick_for_log-domain_calculations
#' @param x a vector of values
#' @return log(sum(exp(x)))
#' @export
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}
