#' Custom rexgauss with default  values
#'
#' @param n Number of sampels to draw, has to be a scalar natural
#' @param mu Mean argument, mu unbound
#' @param sigma Shape parameter, sigma > 0
#' @param beta Shape parameter, beta > 0
#'
#' @return Vector of length n in exgaussian distribution
#' @export
#'
#' @examples hist(rstudent_custom(100, 1, 2, 1))
rexgauss_custom <- function(n, mu = 0, sigma = 1, beta = 1) {
  if(!isNat_len(n)) {
    stop("n has to be a scalar natural number")
  }
  if(isTRUE(any(sigma <= 0))) {
    stop("sigma has to be bigger than 0")
  }
  if(isTRUE(any(beta <= 0))) {
    stop("beta has to be bigger than 0")
  }
  if(!lenEqual(list(mu, sigma, beta), scalars_allowed = TRUE,
               type_check = is.numeric, na_allowed = FALSE)) {
    stop("Type error, or unequal length vectors with n > 1")
  }
  return(brms::rexgaussian(n,  mu, sigma, beta))
}

#' Log-Likelihood of the gumbel distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#'
#' @return Log-Likelihood of gumbel given data in prep
log_lik_exgauss2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  y <- prep$data$Y[i]
  #return(dexgauss2(y, mu, sigma, beta, log = TRUE))
  return(brms::dexgaussian(y, mu, sigma, beta, log = TRUE))
}


#' posterior_predict for the gumbel distribution
#'
#' @param i BRMS indices
#' @param prep BRMS data
#' @param ... catchall argument
#'
#' @return Draws from the Posterior Predictive Distribution
posterior_predict_exgauss2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  #return(rexgauss2(prep$ndraws, mu, sigma, beta))
  return(brms::rexgaussian(prep$ndraws, mu, sigma, beta))
}

#' posterior_epred for the gumbel distribution
#'
#' @param prep BRMS data
#'
#' @return Expected Values of the Posterior Predictive Distribution
posterior_epred_exgauss2 <- function(prep) {
  return(brms::get_dpar(prep, "mu"))
}


#' Exgauss brms custom family, custom implementation vs Stan
#'
#' @param link Link function for function
#' @param link_sigma Link function for sigma argument
#' @param link_beta Link function for beta argument
#'
#' @return BRMS Gumbel distribution family
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = brms::rexgaussian(n = 1000, mu = a + 2, sigma = 2, beta = 2))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'   family = exgauss2(), stanvars = exgauss2()$stanvars,
#'   refresh = 0)
#' plot(fit)
exgauss2 <- function(link = "identity", link_sigma = "log", link_beta = "log") {
  family <- brms::custom_family(
    "exgauss2",
    dpars = c("mu", "sigma", "beta"),
    links = c(link, link_sigma, link_beta),
    lb = c(NA, 0, 0),
    ub = c(NA, NA, NA),
    type = "real",
    log_lik = log_lik_exgauss2,
    posterior_predict = posterior_predict_exgauss2,
    posterior_epred = posterior_epred_exgauss2
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real exgauss2_lpdf(real y, real mu, real sigma, real beta) {
        real y0 = mu - beta;
        real z = y0 + (sigma^2)/beta - y;
        //  lpdf <- log(0.5) - log(beta) + (0.5/beta)*(y0+z-x) + log(2) + log1p(-pnorm(z/sigma))
        return log(0.5)-log(beta)+(0.5/beta)*(y0+z-y)+log(2)+log1m(normal_cdf(z/sigma, 0, 1));
      }

      real exgauss2_rng(real mu, real sigma, real beta) {
        real y0 = mu - beta;
        return normal_rng(y0, sigma) + exponential_rng(1/beta);
      }",
    block = "functions"
  )
  return(family)
}

#' dexgauss2 <- function(x, mu, sigma, beta, log = FALSE) {
#'   # check the arguments
#'   if (isTRUE(any(sigma <= 0))) {
#'     stop("exgauss2 is only defined for sigma > 0")
#'   }
#'   if (isTRUE(any(beta <= 0))) {
#'     stop("exgauss2 is only defined for beta > 0")
#'   }
#'   # Maybe overkill?
#'   if (!lenEqual(list_of_vectors = list(x, mu, sigma), scalars_allowed = TRUE, type_check = is.numeric)) {
#'     stop("exgauss2 argument vectors could not be matched. May be due to wrong type,
#'          or different lengths. Note: len=1 is always allowed, even if the other vectors are len!=1.")
#'   }
#'   if (!isLogic_len(log)) {
#'     stop("the log argument of a density has to be a scalar boolean")
#'   }
#'
#'   # I think, erfc(x) = -2pnorm(x*sqrt(2))
#'   y0 <- mu - beta
#'   z <- y0 + (sigma^2)/beta - x
#'   # lpdf <- log(0.5) - log(beta) + (0.5/beta)*(y0+z-x) - log(2*pnorm(z/sigma))
#'   # lpdf <- log(0.5) - log(beta) + (0.5/beta)*(y0+z-x) + log1p(-erf(z/(sqrt(2)*sigma)))
#'   lpdf <- log(0.5) - log(beta) + (0.5/beta)*(y0+z-x) + log(2) + log1p(-pnorm(z/sigma))
#'   # lpdf <- log(0.5) - log(beta) + (0.5/beta)*(2*(mu-beta) + (sigma*sigma)/beta - 2*x) + log(2) + log1p(-pnorm(mu-beta + (sigma*sigma)/beta - x))
#'
#'   #return either the log or the pdf itself, given the log-value
#'   if (log) {
#'     return(lpdf)
#'   } else {
#'     return(exp(lpdf))
#'   }
#' }
#'
#' rexgauss2 <- function(n, mu = 0, sigma = 1, beta = 1) {
#'   if(!isNat_len(n)) {
#'     stop("n has to be a scalar natural number")
#'   }
#'   if(isTRUE(any(sigma <= 0))) {
#'     stop("sigma has to be bigger than 0")
#'   }
#'   if(isTRUE(any(beta <= 0))) {
#'     stop("beta has to be bigger than 0")
#'   }
#'   if(!lenEqual(list(mu, sigma, beta), scalars_allowed = TRUE,
#'                type_check = is.numeric, na_allowed = FALSE)) {
#'     stop("Type error, or unequal length vectors with n > 1")
#'   }
#'   y0 <- mu - beta
#'   return(rnorm(n, mean = y0, sd = sigma^2) + rexp(n, rate = 1 / beta))
#' }
