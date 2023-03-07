#library(brms)


#pdf
dshifted_inv_gaussian <- function(x, mu = 1, shape = 1, shift = 0, log = FALSE) {
    brms::dinv_gaussian(x-shift, mu, shape, log)
}

#rng
rshifted_inv_gaussian <- function(n, mu = 1, shape = 1, shift = 0) {
    brms::rinv_gaussian(n, mu, shape) + shift
}

posterior_epred_shifted_inverse.gaussian <- function(prep) {
    with(prep$dpars, mu + ndt)
}

posterior_predict_shifted_inverse.gaussian <- function(i, prep, ...) {
    n <- prep$ndraws
    mu <- brms::get_dpar(prep, "mu", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    brms::rshifted_inv_gaussian(n, mu, shape, ndt)
}

log_lik_shifted_inverse.gaussian <- function(i, prep) {
    mu <- brms::get_dpar(prep, "mu", i = i)
    shape <- brms::get_dpar(prep, "shape", i = i)
    ndt <- brms::get_dpar(prep, "ndt", i = i)
    y <- prep$data$Y[i]
    dshifted_inv_gaussian(y, mu, shape, ndt, log = TRUE)
}

shifted_inverse.gaussian <- function(link = "1/mu^2", link_shape = "log", link_ndt = "log"){
    custom_family(
        "shifted_inv_gaussian",
        dpars = c("mu", "shape", "ndt"),
        links = c(link, link_shape, link_ndt),
        lb = c(0, 0, 0),
        ub = c(NA, NA, NA),
        type = "real",
        log_lik = log_lik_shifted_inverse.gaussian,
        posterior_predict = posterior_predict_shifted_inverse.gaussian,
        posterior_epred = posterior_epred_shifted_inverse.gaussian
    )
}

# additionally required Stan code
# stan_shifted_inverse.gaussian <- "
# /* inverse Gaussian log-PDF for a single response
#    * Args:
#    *   y: the response value
#    *   mu: positive mean parameter
#    *   shape: positive shape parameter
#    * Returns:
#    *   a scalar to be added to the log posterior
#    */
#   real inv_gaussian_lpdf(real y, real mu, real shape) {
#     return 0.5 * log(shape / (2 * pi())) - 1.5 * log(y)
#            - 0.5 * shape * square((y - mu) / (mu * sqrt(y)));
#   }
#
# real shifted_inv_gaussian_lpdf(real y, real mu, real shape, real ndt) {
#     return inv_gaussian_lpdf(y - ndt| mu, shape);
# }
# "
#
# stanvars_shifted_inverse.gaussian <- stanvar(scode = stan_shifted_inverse.gaussian, block = "functions")


# ## example
# #simulate some data
# muA = 2
# muB = 3
# shapeA = 4
# shapeB = 5
# ndt = 0.5
# N_per_cond = 1000
# data = data.frame( cond = rep( c('A','B') , each = N_per_cond)
#                    , RT = c( rshifted_inv_gaussian(N_per_cond, muA, shapeA, ndt), rshifted_inv_gaussian(N_per_cond, muB, shapeB, ndt) )
#                   )
#
# library(ggplot2)
# ggplot(data, aes(x=RT, color=cond)) + geom_density()
#
#
# #fit with default link for mu, ie 1/mu^2
# fit1 = brm( bf(RT ~ cond) + shifted_inverse.gaussian()
#            , data
#            , stanvars = stanvars_shifted_inverse.gaussian
#            , cores = 4
# )
# summary(fit1)
# mcmc_plot(fit1)
# pp_check(fit1)
# conditional_effects(fit1)
# conditional_effects(fit1, method = "posterior_predict")
#
# #try with log link for mu, and specifying formulas for shape and ndt parameters
# fit2 = brm( bf(RT ~ cond, shape ~ cond, ndt ~ cond) + shifted_inverse.gaussian("log")
#             , data
#             , stanvars = stanvars_shifted_inverse.gaussian
#             , prior = prior(normal(0, 1) , class = "b")
#             , cores = 4
# )
# summary(fit2)
# mcmc_plot(fit2)
# pp_check(fit2)
# conditional_effects(fit2)
# conditional_effects(fit2, method = "posterior_predict")

