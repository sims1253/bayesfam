#' Lookup function for brms families via string identifier
#'
#' @param family String identifier of the family.
#' @param link Link to be passed to the family function.
#'
#' @return A brmsfamily object matching the string identifier and using the link
#' @export
#'
#' @examples
#' brms_family_lookup("weibull", "softplus")
brms_family_lookup <- function(family, link = NULL) {
  switch(family,
    "beta" = brms::brmsfamily("beta", link = link),
    "kumaraswamy" = bayesfam::kumaraswamy(link = link),
    "logitnormal" = bayesfam::logitnormal(link = link),
    "cauchitnormal" = bayesfam::cauchitnormal(link = link),
    "cloglognormal" = bayesfam::cloglognormal(link = link),
    "simplex" = bayesfam::simplex(link = link),
    "gaussian" = brms::brmsfamily("gaussian", link = link),
    "gamma" = brms::brmsfamily("gamma", link = link),
    "weibull" = brms::brmsfamily("weibull", link = link),
    "lognormal" = brms::brmsfamily("lognormal", link = link),
    "softplusnormal" = bayesfam::softplusnormal(link = link),
    "lomax" = bayesfam::lomax(link = link),
    "frechet" = brms::brmsfamily("frechet", link = link),
    "inverse.gaussian" = brms::brmsfamily("inverse.gaussian", link = link),
    "betaprime" = bayesfam::betaprime(link = link),
    "gompertz" = bayesfam::gompertz(link = link)
  )
}

#' Lookup function for RNGs via string identifier
#'
#' @param family String identifier of the likelihood family to get an RNG for.
#'
#' @return The RNG function.
#' @export
#'
#' @examples
#' rng_lookup("gamma")
#' do.call(rng_lookup("gaussian"), list(n = 100, mean = 0, sd = 1))
rng_lookup <- function(family) {
  switch(family,
    "beta" = bayesfam::rbeta_mean,
    "kumaraswamy" = bayesfam::rkumaraswamy,
    "logitnormal" = bayesfam::rlogitnormal,
    "cauchitnormal" = bayesfam::rcauchitnormal,
    "cloglognormal" = bayesfam::rcloglognormal,
    "simplex" = bayesfam::rsimplex,
    "gaussian" = rnorm,
    "gamma" = bayesfam::rgamma_mean,
    "weibull" = bayesfam::rweibull_median,
    "lognormal" = bayesfam::rlognormal,
    "softplusnormal" = bayesfam::rsoftplusnormal,
    "lomax" = bayesfam::rlomax,
    "frechet" = bayesfam::rfrechet_median,
    "inverse.gaussian" = brms::rinv_gaussian,
    "betaprime" = bayesfam::rbetaprime,
    "gompertz" = bayesfam::rgompertz
  )
}


#' Lookup function for link and repsonse functions via string identifier.
#'
#' If a transformed normal likelihood is passed, the respective built-in link will
#' be returned instead of the identity link that would commonly be used with
#' transformed normal likelihood families.
#'
#' @param link String identifier for the link function of interest.
#' @param family If a transformed normal family is passed, returns the
#'               respective link instead of \code{link}
#' @param inv True to return the response function instead of the link function.
#'
#' @return The respective link function.
#' @export
#'
#' @examples
#'
#' link_lookup("log", "gaussian", FALSE)
#'
#' link_lookup("identiy", "logitnormal", FALSE)
link_lookup <- function(link, family = NULL, inv = FALSE) {
  if (inv) {
    if (!is.null(family)) {
      switch(family,
        "logitnormal" = return(bayesfam::inv_logit),
        "cauchitnormal" = return(bayesfam::inv_cauchit),
        "cloglognormal" = return(bayesfam::inv_cloglog),
        "lognormal" = return(exp),
        "softplusnormal" = return(bayesfam::inv_softplus)
      )
    }
    switch(link,
      "logit" = bayesfam::inv_logit,
      "cauchit" = bayesfam::inv_cauchit,
      "cloglog" = bayesfam::inv_cloglog,
      "identity" = identity,
      "log" = exp,
      "softplus" = bayesfam::inv_softplus
    )
  } else {
    if (!is.null(family)) {
      switch(family,
        "logitnormal" = return(bayesfam::logit),
        "cauchitnormal" = return(bayesfam::cauchit),
        "cloglognormal" = return(bayesfam::cloglog),
        "lognormal" = return(log),
        "softplusnormal" = return(bayesfam::softplus)
      )
    }
    switch(link,
      "logit" = bayesfam::logit,
      "cauchit" = bayesfam::cauchit,
      "cloglog" = bayesfam::cloglog,
      "identity" = identity,
      "log" = log,
      "softplus" = bayesfam::softplus
    )
  }
}

#' Lookup function for the names of the auxiliary parameters of a likelihood
#'
#' @param family
#'
#' @return
#' @export
#'
#' @examples
aux_family_parameters_lookup <- function(family) {
  brms_family_lookup(family)$dpars[2:length(brms_family_lookup(family)$dpars)]
}

#' Lookup for limits of family auxiliary parameters.
#'
#' @param family The identifier string of a family.
#'
#' @return List containing lower and upper bounds for the auxiliary parameter.
#' @export
#'
#' @examples
aux_limits_lookup <- function(family) {
  switch(family,
    "beta" = list(lb = 0, ub = Inf),
    "kumaraswamy" = list(lb = 0, ub = Inf),
    "logitnormal" = list(lb = 0, ub = Inf),
    "cauchitnormal" = list(lb = 0, ub = Inf),
    "cloglognormal" = list(lb = 0, ub = Inf),
    "simplex" = list(lb = -Inf, ub = Inf),
    "gaussian" = list(lb = 0, ub = Inf),
    "gamma" = list(lb = 0, ub = Inf),
    "weibull" = list(lb = 0, ub = Inf),
    "lognormal" = list(lb = 0, ub = Inf),
    "softplusnormal" = list(lb = 0, ub = Inf),
    "lomax" = list(lb = 1, ub = Inf),
    "frechet" = list(lb = 1, ub = Inf),
    "inverse.gaussian" = list(lb = 0, ub = Inf),
    "betaprime" = list(lb = 0, ub = Inf),
    "gompertz" = list(lb = 0, ub = Inf)
  )
}
