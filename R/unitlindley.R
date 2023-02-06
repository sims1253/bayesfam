# suggested in doi:10.1080/02664763.2018.1511774 for data on the unit interval
# https://www.tandfonline.com/doi/full/10.1080/02664763.2018.1511774

#library(brms)
#library(lamW) # for the quantile function in posterior_predict

# helper functions for post-processing of the family
log_lik_unitlindley <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  y <- prep$data$Y[i]
  2*log(1-mu)-log(mu)-3*log(1-y)-y*(1-mu)/(mu*(1-y))
}

posterior_epred_unitlindley <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

posterior_predict_unitlindley <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  lambert <- lamW::lambertWm1((runif(1, 0, 1)-1)/mu*exp(-1/mu))
  return((1/mu+lambert)/(1+lambert))
}

dunit_lindley <- function(x, mu, log=FALSE) {
  if(isTRUE(any(x <= 0 | x >= 1))) {
    stop("The x argument has to be in the interval (0, 1)")
  }
  if(isTRUE(any(mu <= 0 | mu >= 1))) {
    stop("The mu argument has to be in the interval (0, 1)")
  }
  if(!isLogic_len(log)) {
    stop("The log argument has to be a single boolean value")
  }

  lpdf <- 2*log(1-mu)-log(mu)-3*log(1-x)-x*(1-mu)/(mu*(1-x))
  if(log) {
    return(lpdf)
  }
  else {
    return(exp(lpdf))
  }
}

# theta = 1/mu - 1

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

runit_lindley <- function(n, mu) {
  if(!isNat_len(n)) {
    stop("n has to be a single natural number")
  }
  return(qunit_lindley(runif(n), mu))
}

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
        return 2*log(1-mu)-log(mu)-3*log(1-y)-y*(1-mu)/(mu*(1-y));
      }
    ",
    block = "functions"
  )
}
