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
  dbeta_binomial2(y, mu, phi, trials)
}

posterior_predict_beta_binomial2 <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  trials <- prep$data$vint1[i]
  rbeta_binomial2(mu, phi, trials)
}

posterior_epred_beta_binomial2 <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  trials <- prep$data$vint1
  trials <- matrix(trials, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE)
  mu * trials
}

# definition of the custom family
#' Title
#'
#' @param link
#' @param link_phi
#'
#' @return
#' @export
#'
#' @examples
beta_binomial2 <- function(link = "logit", link_phi = "log") {
  family <- brms::custom_family(
    "beta_binomial2",
    dpars = c("mu", "phi"),
    links = c(link, link_phi),
    lb = c(0, 0),
    ub = c(1, NA),
    type = "int",
    vars = "trials[n]",
    log_lik = log_lik_beta_binomial2,
    posterior_predict = posterior_predict_beta_binomial2,
    posterior_epred = posterior_epred_beta_binomial2
  )
  family$stanvars <- brms::stanvar(
    scode = "
      real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
        return(beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi));
      }

      int beta_binomial2_rng(real mu, real phi, int T) {
        return(beta_binomial_rng(T, mu * phi, (1 - mu) * phi));
      }",
    block = "functions"
  )
}

# family$stanvars <- brms::stanvar(
#   scode = "
#       real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
#         return(beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi));
#       }
#       int beta_binomial2_rng(real mu, real phi, int T) {
#         return(beta_binomial_rng(T, mu * phi, (1 - mu) * phi));
#       }",
#   block = "functions"
# )

#' Title
#'
#' @param x
#' @param mu
#' @param phi
#' @param t
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dbeta_binomial2 <- function(x, mu, phi, t, log = FALSE) {
  # check the arguments
  if (!isNat_len(x, len=length(x))) {
    stop("beta binomial2 is only defined for x > 0")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("beta binomial2 is only defined for phi > 0")
  }
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop("beta binomial2 is only defined for mu e [0, 1]")
  }
  if (!isInt_len(t, len=length(t))) {
    stop("Only natural numbers are allowed for the t argument in beta binomial2")
  }
  if (isFALSE(isLogic_len(log))) {
    stop("The log argument has to be a boolean of len 1 in beta binomial2")
  }
  if (isTRUE(any(x > t))) {
    stop("No x may be bigger than the t argument in beta binomial2")
  }

  # calculate the usual alpha and beta for the beta function
  a <- mu*phi
  b <- (1-mu)*phi
  # calculate stability optimized log of the pdf
  lpdf <- log_bin_coeff(t,x) + log_beta(x + a, t - x + b) - log_beta(a, b)
  #lpdf <- log(choose(t,x)) + log(beta(x + a, t - x + b)) - log(beta(a, b))

  # return the corresponding value
  if(log) {
    return(lpdf)
  } else {
    return(exp(lpdf))
  }
}

log_beta <- function(a, b) {
  # Is this actually better than log(beta)?
  return(log(gamma(a)) + log(gamma(b)) - log(gamma(a + b)))
}

log_bin_coeff <- function(n, k) {
  # Is this actually better, than log(choose)?
  # https://en.wikipedia.org/wiki/Binomial_coefficient
  # log of bin coefficient, to be more stable, with sum, rather than mult
  #return(sum_nat_numbers(n) - sum_nat_numbers(k) - sum_nat_numbers(n-k))
  return(log_fact(n) - log_fact(k) - log_fact(n-k))
}

log_fact <- function(n) {
  # Using Pochhammer with x=1
  # log(n!) = sum log(n)
  # -> https://www.wolframalpha.com/input?i=sum+log%28n%29
  return(log(gamma(1 + n)))
}

# pochhammer <- function(x, n) {
#   # https://mathworld.wolfram.com/PochhammerSymbol.html
# }

# pochhammer_x1 <- function(n) {
#   # for x=1, the formula gets simpler
#   # https://mathworld.wolfram.com/PochhammerSymbol.html
#   return(gamma(1+n))
# }

# beta2 <- function(mu, phi) {
#   return(beta(mu*phi, (1-mu)*phi))
# }

#' Title
#'
#' @param n
#' @param mu
#' @param phi
#' @param t
#'
#' @return
#' @export
#'
#' @examples
rbeta_binomial2 <- function(n, mu, phi, t) {
  # check the n argument, rest is checked by quantile function
  if (!isNat_len(n)) {
    stop("In rbeta_binomial2 the n argument has to be a natural scalar")
  }
  if (isTRUE(any(phi <= 0))) {
    stop("beta binomial2 is only defined for phi > 0")
  }
  if (isTRUE(any(mu < 0 | mu > 1))) {
    stop("beta binomial2 is only defined for mu e [0, 1]")
  }
  if (!isInt_len(t, len=length(t))) {
    stop("Only natural numbers are allowed for the t argument in beta binomial2")
  }

  # calculate beta_binom, by calculating beta2 and using the result for binom
  p <- rbeta(n, mu*phi, (1-mu)*phi)
  return(rbinom(n, t, p))
}

# additionally required Stan code
# stan_beta_binomial2 <- "
#   real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
#     return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
#   }
#   int beta_binomial2_rng(real mu, real phi, int T) {
#     return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
#   }
# "
