
# helper functions for post-processing of the family
log_lik_lognormal_natural <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)

  # if(NCOL(brms::get_dpar(prep, "sigma"))==1){
  #   sigma <- brms::get_dpar(prep, "sigma")
  # }else{
  #   # is this really necassary?
  #   sigma <- brms::get_dpar(prep, "sigma", i=i)
  # }

  y <- prep$data$Y[i]
  #common_term = log1p(sigma^2/mu^2)
  #Vectorize(dlnorm)(y, log(mu)-common_term/2, sqrt(common_term), log = TRUE)
  return(dlognormal_natural(y, mu, sigma, log=TRUE))
}


posterior_predict_lognormal_natural <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)

  # if(NCOL(brms::get_dpar(prep, "sigma"))==1){
  #   sigma <- brms::get_dpar(prep, "sigma")
  # }else{
  #   # is this really necassary?
  #   sigma <- brms::get_dpar(prep, "sigma", i=i)
  # }   ## [, i] if sigma is modelled, without otherwise

  #common_term = log(1+sigma^2/mu^2)
  #rlnorm(n, log(mu)-common_term/2, sqrt(common_term))
  return(rlognormal_natural(prep$ndraws, mu, sigma))
}

posterior_epred_lognormal_natural <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  return(mu)
}


#' Title
#'
#' @param x Value space, x > 0
#' @param mu Mean parameter, mu > 0
#' @param sigma Shape parameter, sigma unbound
#' @param log boolean value to return log-pdf, default FALSE
#'
#' @return
#' @export
#'
#' @examples
dlognormal_natural <- function(x, mu, sigma, log = FALSE) {
  if(isTRUE(any(mu <= 0))) {
    stop("Mu has to be > 0")
  }
  common_term <- log1p(sigma^2/mu^2)
  return(dlognormal(x, log(mu)-common_term/2, sqrt(common_term), log))
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
rlognormal_natural <- function(n, mu = 0, sigma = 1) {
  if(isTRUE(any(mu <= 0))) {
    stop("Mu has to be > 0")
  }
  common_term <- log1p(sigma^2/mu^2)
  return(rlognormal(n, log(mu)-common_term/2, sqrt(common_term)))
}

#' Title
#'
#' @param link
#' @param link_sigma
#'
#' @return
#' @export
#'
#' @examples a <- rnorm(n = 1000)
#' data <- list(a = a, y = rlognormal_natural(n = 1000, mu = exp(0.5 * a + 1), sigma = 2))
#' fit <- brms::brm(formula = y ~ 1 + a, data = data,
#'  family = lognormal_natural(), stanvars = lognormal_natural()$stanvars,
#'  refresh = 0)
#' plot(fit)
lognormal_natural <- function(link = "log", link_sigma = "log") {
  family <- brms::custom_family(
    "lognormal_natural",
    dpars = c("mu", "sigma"),
    links = c(link, link_sigma),
    lb = c(0, 0),
    type = "real",
    log_lik = log_lik_lognormal_natural,
    posterior_predict = posterior_predict_lognormal_natural,
    posterior_epred = posterior_epred_lognormal_natural
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real lognormal_natural_lpdf(real y, real mu, real sigma) {
        real common_term = log1p(sigma^2/mu^2);
        return lognormal_lpdf(y | log(mu)-common_term/2,
                                  sqrt(common_term));
      }
      real lognormal_natural_rng(real mu, real sigma) {
        real common_term = log1p(sigma^2/mu^2);
        return lognormal_rng(log(mu)-common_term/2,
                                sqrt(common_term));
      }
    ",
    block = "functions"
  )
  return(family)
}
