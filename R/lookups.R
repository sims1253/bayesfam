# The single source of per-family metadata behind the lookup functions below.
# What can be derived from the family object (auxiliary parameter names,
# bounds of custom families) is derived instead of stored. Entries hold:
# * family: constructor, called without arguments to keep its defaults.
# * rng: RNG in the family's parameterization.
# * lb/ub: auxiliary parameter bounds; only stored for brms built-ins.
# * link: built-in link of transformed normal families (see link_lookup()).
# A function rather than a list so that it can refer to functions defined in
# files collated later.
family_registry <- function() {
  custom <- function(family, rng, link = NULL) {
    list(family = family, rng = rng, link = link)
  }
  builtin <- function(name, rng, lb, ub = Inf, link = NULL) {
    family <- function(link = NULL) brms::brmsfamily(name, link = link)
    list(family = family, rng = rng, lb = lb, ub = ub, link = link)
  }
  list(
    "beta" = builtin("beta", rbeta_mean, lb = 0),
    "betaprime" = custom(betaprime, rbetaprime),
    "cauchitnormal" = custom(cauchitnormal, rcauchitnormal, "cauchit"),
    "cloglognormal" = custom(cloglognormal, rcloglognormal, "cloglog"),
    # brms frechet uses mu = mean, which requires nu > 1
    "frechet" = builtin(
      "frechet",
      function(n, mu, nu) {
        if (isTRUE(any(nu <= 1))) {
          stop("frechet is only defined for nu > 1")
        }
        brms::rfrechet(n, scale = mu / gamma(1 - 1 / nu), shape = nu)
      },
      lb = 1
    ),
    "gamma" = builtin("gamma", rgamma_mean, lb = 0),
    "gaussian" = builtin("gaussian", rnorm, lb = 0),
    "generalized_gamma" = custom(generalized_gamma, rgeneralized_gamma),
    "generalized_normal" = custom(generalized_normal, rgeneralized_normal),
    "gompertz" = custom(gompertz, rgompertz),
    "gumbel_mean" = custom(gumbel_mean, rgumbel_mean),
    "inverse.gaussian" = builtin(
      "inverse.gaussian",
      brms::rinv_gaussian,
      lb = 0
    ),
    "kumaraswamy" = custom(kumaraswamy, rkumaraswamy),
    "logistic" = custom(logistic, rlogistic),
    "logitnormal" = custom(logitnormal, rlogitnormal, "logit"),
    "lognormal" = builtin("lognormal", rlognormal, lb = 0, link = "log"),
    "lognormal_natural" = custom(lognormal_natural, rlognormal_natural),
    "lomax" = custom(lomax, rlomax),
    "shifted_inv_gaussian" = custom(
      shifted_inv_gaussian,
      rshifted_inv_gaussian
    ),
    "shifted_lognormal_uniform" = custom(
      shifted_lognormal_uniform,
      rshifted_lognormal_uniform
    ),
    "simplex" = custom(simplex, rsimplex),
    "softplusnormal" = custom(softplusnormal, rsoftplusnormal, "softplus"),
    "symlognormal" = custom(symlognormal, rsymlognormal, "symlog"),
    "unit_lindley" = custom(unit_lindley, runit_lindley),
    # brms weibull uses mu = mean
    "weibull" = builtin(
      "weibull",
      function(n, mu, shape) rweibull(n, shape, mu / gamma(1 + 1 / shape)),
      lb = 0
    )
  )
}

# Registry entry of `family`, erroring on unknown identifiers.
family_entry <- function(family) {
  registry <- family_registry()
  if (!isSingleString(family) || !family %in% names(registry)) {
    stop(
      "Unknown family '",
      format(family),
      "'. Supported families: ",
      paste(names(registry), collapse = ", "),
      call. = FALSE
    )
  }
  registry[[family]]
}

#' Lookup function for brms families via string identifier
#'
#' @param family String identifier of the family. Unknown identifiers raise an
#'   error.
#' @param link Link to be passed to the family function. If `NULL`, the
#'   family's default link is used.
#'
#' @return A brmsfamily object matching the string identifier and using the link
#' @export
#'
#' @examples
#' brms_family_lookup("weibull", "softplus")
brms_family_lookup <- function(family, link = NULL) {
  constructor <- family_entry(family)$family
  if (is.null(link)) {
    return(constructor())
  }
  constructor(link = link)
}

#' Lookup function for RNGs via string identifier
#'
#' @param family String identifier of the likelihood family to get an RNG for.
#'   Unknown identifiers raise an error.
#'
#' @return The RNG function.
#' @export
#'
#' @examples
#' rng_lookup("gamma")
#' do.call(rng_lookup("gaussian"), list(n = 100, mean = 0, sd = 1))
rng_lookup <- function(family) {
  family_entry(family)$rng
}


#' Lookup function for link and repsonse functions via string identifier.
#'
#' If a transformed normal likelihood is passed, the respective built-in link will
#' be returned instead of the identity link that would commonly be used with
#' transformed normal likelihood families.
#'
#' @param link String identifier for the link function of interest. Unknown
#'   identifiers raise an error.
#' @param family If a transformed normal family is passed, returns the
#'               respective link instead of `link`
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
  if (isSingleString(family)) {
    family_link <- family_registry()[[family]]$link
    if (!is.null(family_link)) {
      link <- family_link
    }
  }
  links <- list(
    "logit" = list(logit, inv_logit),
    "cauchit" = list(cauchit, inv_cauchit),
    "cloglog" = list(cloglog, inv_cloglog),
    "identity" = list(identity, identity),
    "log" = list(log, exp),
    "softplus" = list(softplus, inv_softplus),
    "symlog" = list(symlog, inv_symlog)
  )
  if (!isSingleString(link) || !link %in% names(links)) {
    stop(
      "Unknown link '",
      format(link),
      "'. Supported links: ",
      paste(names(links), collapse = ", "),
      call. = FALSE
    )
  }
  links[[link]][[1 + inv]]
}

#' Lookup function for the names of the auxiliary parameters of a likelihood
#'
#' @param family The identifier string of a family. Unknown identifiers raise
#'   an error.
#'
#' @return Character vector of auxiliary parameter names (all distributional
#'   parameters except `mu`), empty if there are none.
#' @export
#'
#' @examples
#' aux_family_parameters_lookup("beta")
aux_family_parameters_lookup <- function(family) {
  setdiff(brms_family_lookup(family)$dpars, "mu")
}

#' Lookup for limits of family auxiliary parameters.
#'
#' @param family The identifier string of a family. Unknown identifiers raise
#'   an error.
#'
#' @return List containing lower and upper bounds for the auxiliary parameters,
#'   in the order of [aux_family_parameters_lookup()].
#' @export
#'
#' @examples
#' aux_limits_lookup("beta")
aux_limits_lookup <- function(family) {
  entry <- family_entry(family)
  if (!is.null(entry$lb)) {
    return(list(lb = entry$lb, ub = entry$ub))
  }
  family <- entry$family()
  aux <- family$dpars != "mu"
  # unbounded (NA) and data-dependent (e.g. "min_Y") bounds become infinite
  bound <- function(x, default) {
    x <- suppressWarnings(as.numeric(unlist(x[aux], use.names = FALSE)))
    replace(x, is.na(x), default)
  }
  list(lb = bound(family$lb, -Inf), ub = bound(family$ub, Inf))
}
