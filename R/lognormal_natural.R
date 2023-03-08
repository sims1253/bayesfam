
# helper functions for post-processing of the family
log_lik_lognormal_natural <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)

  if(NCOL(brms::get_dpar(prep, "sigma"))==1){
    sigma <- brms::get_dpar(prep, "sigma")
  }else{
    # is this really necassary?
    sigma <- brms::get_dpar(prep, "sigma", i=i)
  }

  y <- prep$data$Y[i]
  common_term = log(1+sigma^2/mu^2)
  Vectorize(dlnorm)(y, log(mu)-common_term/2, sqrt(common_term), log = TRUE)
}


posterior_predict_lognormal_natural <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)

  if(NCOL(brms::get_dpar(prep, "sigma"))==1){
    sigma <- brms::get_dpar(prep, "sigma")
  }else{
    # is this really necassary?
    sigma <- brms::get_dpar(prep, "sigma", i=i)
  }   ## [, i] if sigma is modelled, without otherwise

  common_term = log(1+sigma^2/mu^2)
  rlnorm(n, log(mu)-common_term/2, sqrt(common_term))
}

posterior_epred_lognormal_natural <- function(prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  return(mu)
}


dlognormal_natural <- function(x, mu, sigma, log = FALSE) {
  common_term <- log(1+sigma^2/mu^2)
  return(dlognormal(y, log(mu)-common_term/2, sqrt(common_term), log))
}

rlognormal_natural <- function(n, mu, sigma) {
  common_term <- log(1+sigma^2/mu^2)
  return(rlognormal(n, log(mu)-common_term/2, sqrt(common_term)))
}

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
        real common_term = log(1+sigma^2/mu^2);
        return lognormal_lpdf(y | log(mu)-common_term/2,
                                  sqrt(common_term));
      }
      real lognormal_natural_rng(real mu, real sigma) {
        real common_term = log(1+sigma^2/mu^2);
        return lognormal_rng(log(mu)-common_term/2,
                                sqrt(common_term));
      }
    ",
    block = "functions"
  )
  return(family)
}
