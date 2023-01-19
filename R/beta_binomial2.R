# This is the beta-binomial distribution in mean-precision parameterization
# Details are provided in
# https://paul-buerkner.github.io/brms/articles/brms_customfamilies.html

library(brms)

# helper functions for post-processing of the family
log_lik_beta_binomial2 <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  y <- prep$data$Y[i]
  beta_binomial2_lpmf(y, mu, phi, trials)
}

posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  beta_binomial2_rng(mu, phi, trials)
}

posterior_epred_beta_binomial2 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  trials <- prep$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

# definition of the custom family
beta_binomial2 <- function(link = "logit", link_phi = "log") {
  custom_family(
    "beta_binomial2",
    dpars = c("mu", "phi"),
    links = c(link, link_phi),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "int",
    vars = "trials[n]",
    log_lik = log_lik_beta_binomial2,
    posterior_predict = posterior_predict_beta_binomial2,
    posterior_epred = posterior_epred__beta_binomial2
  )
}

dbeta_binomial2 <- function(x, mu, phi, t, log = FALSE) {
  # check the arguments
  if (!isNat_len(x, len=length(x))) {
    stop("beta prime is only defined for x > 0")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("beta prime is only defined for phi > 0")
  }
  if (isTRUE(any(mu <= 0))) {
    stop("beta prime is only defined for mu > 0")
  }
  if (!isInt_len(t, len=length(t))) {
    stop("Only natural numbers are allowed for the t argument")
  }
  if (isFALSE(isLogic_len(log))) {
    stop("The log argument has to be a boolean of len 1")
  }

  #lpdf <- beta_binomial2_lpmf(x, mu, phi, t)
  a <- mu*phi
  b <- (1-mu)*phi
  lpdf <- log_bin_coeff(t, x) + log(beta(x + a, t - x + b)) - log(beta(a, b))
  if(log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

sum_nat_numbers <- function(n) {
  # small helpers, n is already checked by caller here, so no further checks necassary
  # but I guess, n > 0!
  return(n*(n+1)/2)
}

log_bin_coeff <- function(n, k) {
  # https://en.wikipedia.org/wiki/Binomial_coefficient
  # log of bin coefficient, to be more stable, with sum, rather than mult
  return(sum_nat_numbers(n) - sum_nat_numbers(k) - sum_nat_numbers(n-k))
}

beta2 <- function(mu, phi) {
  return(beta(mu*phi, (1-mu)*phi))
}

rbeta_binomial2 <- function(x, mu, phi, t, log = FALSE) {
  # check the arguments
  if (isTRUE(any(x <= 0))) {
    stop("beta prime is only defined for x > 0")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("beta prime is only defined for phi > 0")
  }
  if (isTRUE(any(mu <= 0))) {
    stop("beta prime is only defined for mu > 0")
  }
  if (!isNat_len(t, len=length(t))) {
    stop("Only natural numbers are allowed for the t argument")
  }
  if (isFALSE(isLogic_len(log))) {
    stop("The log argument has to be a boolean of len 1")
  }

  return(beta_binomial2_rng(mu, phi, t))
}

# additionally required Stan code
stan_beta_binomial2 <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
