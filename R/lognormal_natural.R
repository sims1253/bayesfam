# log-normal distribution in natural scale (mean-sd) parameterization

library(brms)

# helper functions for post-processing of the family
log_lik_lognormal_natural <- function(i, prep) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$sigma)==1){sigma <- prep$dpars$sigma}else
  {sigma <- prep$dpars$sigma[, i]}   ## [, i] if sigma is modelled, without otherwise
  y <- prep$data$Y[i]
  common_term = log(1+sigma^2/mu^2)
  Vectorize(dlnorm)(y, log(mu)-common_term/2, sqrt(common_term), log = TRUE)
}


posterior_predict_lognormal_natural <- function(i, prep, ...) {
  mu <- prep$dpars$mu[, i]
  if(NCOL(prep$dpars$sigma)==1){sigma <- prep$dpars$sigma}else
  {sigma <- prep$dpars$sigma[, i]}   ## [, i] if sigma is modelled, without otherwise
  common_term = log(1+sigma^2/mu^2)
  rlnorm(n, log(mu)-common_term/2, sqrt(common_term))
}

posterior_epred_lognormal_natural <- function(prep) {
  mu <- prep$dpars$mu
  return(mu)
}

# definition of the custom family
custom_family(name = "lognormal_natural",
              dpars = c("mu", "sigma"),
              links = c("log", "log"),
              lb = c(0, 0),
              type = "real") ->
  lognormal_natural

# additionally required Stan code
stan_lognormal_natural <- "
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
"


# # example
# stanvars <- stanvar(scode = stan_lognormal_natural, block = "functions")
#
# brm(data = mtcars, formula = bf(mpg ~ wt),
#     family = lognormal_natural, stanvars = stanvars) -> model_test1
#
# brm(data = mtcars, formula = bf(mpg ~ wt, sigma ~ wt),
#     family = lognormal_natural, stanvars = stanvars) -> model_test2
#
# model_test1 |> conditional_effects()
# model_test2 |> conditional_effects()
#
# model_test1 |> conditional_effects(method = "posterior_predict")
# model_test2 |> conditional_effects(method = "posterior_predict")
#
# loo(model_test1, model_test2)
