#' Title
#'
#' @param family
#' @param link
#'
#' @return
#' @export
#'
#' @examples
brms_family_lookup <- function(family, link = NULL) {
  switch(family,
    "beta" = brms::brmsfamily("beta", link = link),
    "kumaraswamy" = kumaraswamy(link = link),
    "logitnormal" = logitnormal(link = link),
    "cauchitnormal" = cauchitnormal(link = link),
    "cloglognormal" = cloglognormal(link = link),
    "simplex" = simplex(link = link),
    "gaussian" = brms::brmsfamily("gaussian", link = link),
    "gamma" = brms::brmsfamily("gamma", link = link),
    "weibull" = brms::brmsfamily("weibull", link = link),
    "lognormal" = brms::brmsfamily("lognormal", link = link),
    "softplusnormal" = softplusnormal(link = link),
    "lomax" = lomax(link = link),
    "frechet" = brms::brmsfamily("frechet", link = link),
    "inverse.gaussian" = brms::brmsfamily("inverse.gaussian", link = link),
    "betaprime" = betaprime(link = link),
    "gompertz" = gompertz(link = link),
    "inversegaussian_custom" = inversegaussian_custom(link = link)
  )
}

#' Title
#'
#' @param family
#' @param link
#'
#' @return
#' @export
#'
#' @examples
rng_lookup <- function(family) {
  switch(family,
    "beta" = rbeta_custom,
    "kumaraswamy" = rkumaraswamy,
    "logitnormal" = rlogitnormal,
    "cauchitnormal" = rcauchitnormal,
    "cloglognormal" = rcloglognormal,
    "simplex" = rsimplex,
    "gaussian" = rnorm,
    "gamma" = rgamma_custom,
    "weibull" = rweibull_custom,
    "lognormal" = rlognormal,
    "softplusnormal" = rsoftplusnormal,
    "lomax" = rlomax,
    "frechet" = rfrechet_custom,
    "inverse.gaussian" = brms::rinv_gaussian,
    "betaprime" = rbetaprime,
    "gompertz" = rgompertz,
    "inversegaussian_custom" = rinversegaussian_custom
  )
}


#' Title
#'
#' @param link
#'
#' @return
#' @export
#'
#' @examples
inv_link_lookup <- function(link) {
  switch(link,
    "logit" = inv_logit,
    "cauchit" = inv_cauchit,
    "cloglog" = inv_cloglog,
    "identity" = identity,
    "log" = exp,
    "softplus" = inv_softplus
  )
}

#' Title
#'
#' @param family
#'
#' @return
#' @export
#'
#' @examples
second_family_parameter_lookup <- function(family) {
  brms_family_lookup(family)$dpars[2]
}

#' Title
#'
#' @param family
#'
#' @return
#' @export
#'
#' @examples
prior_lookup <- function(family) {
  switch(family,
    "frechet" = c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = "nu", lb = 1.00001)
    ),
    c(
      brms::set_prior("", class = "Intercept"),
      brms::set_prior("", class = second_family_parameter_lookup(family))
    )
  )
}
