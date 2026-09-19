# Single family registry powering all lookup helpers.
#
# Each entry provides:
# * constructor: builds the brms family; takes a single optional `link`
#   argument for the mu link (auxiliary links keep their constructor
#   defaults). Called without arguments, constructor defaults apply.
# * rng: the RNG matching the family parameterization.
# * aux: names of the auxiliary distributional parameters (dpars minus mu).
# * aux_lb / aux_ub: numeric bounds of the auxiliary parameters.
# * link: optional override of the natural link for transformed normal
#   families used by link_lookup().
#
# Registry state notes:
# * weibull/frechet map to the median parameterized helpers (#38).
# * simplex has sigma > 0 with a log link (final state of #27).
family_registry <- list(
  "beta" = list(
    constructor = function(link = NULL) brms::brmsfamily("beta", link = link),
    rng = rbeta_mean,
    aux = "phi",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "betaprime" = list(
    constructor = betaprime,
    rng = rbetaprime,
    aux = "phi",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "cauchitnormal" = list(
    constructor = cauchitnormal,
    rng = rcauchitnormal,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf,
    link = "cauchit"
  ),
  "cloglognormal" = list(
    constructor = cloglognormal,
    rng = rcloglognormal,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf,
    link = "cloglog"
  ),
  "frechet" = list(
    constructor = function(link = NULL) brms::brmsfamily("frechet", link = link),
    rng = rfrechet_median,
    aux = "nu",
    aux_lb = 1,
    aux_ub = Inf
  ),
  "gamma" = list(
    constructor = function(link = NULL) brms::brmsfamily("gamma", link = link),
    rng = rgamma_mean,
    aux = "shape",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "gaussian" = list(
    constructor = function(link = NULL) brms::brmsfamily("gaussian", link = link),
    rng = rnorm,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "generalized_gamma" = list(
    constructor = generalized_gamma,
    rng = rgeneralized_gamma,
    aux = c("sigma", "Q"),
    aux_lb = c(0, -Inf),
    aux_ub = c(Inf, Inf)
  ),
  "generalized_normal" = list(
    constructor = generalized_normal,
    rng = rgeneralized_normal,
    aux = c("sigma", "beta"),
    aux_lb = c(0, 0),
    aux_ub = c(Inf, Inf)
  ),
  "gompertz" = list(
    constructor = gompertz,
    rng = rgompertz,
    aux = "beta",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "gumbel_mean" = list(
    constructor = gumbel_mean,
    rng = rgumbel_mean,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "inverse.gaussian" = list(
    constructor = function(link = NULL) brms::brmsfamily("inverse.gaussian", link = link),
    rng = brms::rinv_gaussian,
    aux = "shape",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "inverse_burr" = list(
    constructor = inverse_burr,
    rng = rinverse_burr,
    aux = c("tau", "gamma"),
    aux_lb = c(0, 0),
    aux_ub = c(Inf, Inf)
  ),
  "kumaraswamy" = list(
    constructor = kumaraswamy,
    rng = rkumaraswamy,
    aux = "p",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "logistic" = list(
    constructor = logistic,
    rng = rlogistic,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "logitnormal" = list(
    constructor = logitnormal,
    rng = rlogitnormal,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf,
    link = "logit"
  ),
  "lognormal" = list(
    constructor = function(link = NULL) brms::brmsfamily("lognormal", link = link),
    rng = rlognormal,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf,
    link = "log"
  ),
  "lognormal_natural" = list(
    constructor = lognormal_natural,
    rng = rlognormal_natural,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "lomax" = list(
    constructor = lomax,
    rng = rlomax,
    aux = "alpha",
    aux_lb = 1,
    aux_ub = Inf
  ),
  "shifted_inv_gaussian" = list(
    constructor = shifted_inv_gaussian,
    rng = rshifted_inv_gaussian,
    aux = c("shape", "ndt"),
    # the true upper bound of ndt is the minimum of the data (min_Y)
    aux_lb = c(0, 0),
    aux_ub = c(Inf, Inf)
  ),
  "shifted_lognormal_uniform" = list(
    constructor = shifted_lognormal_uniform,
    rng = rshifted_lognormal_uniform,
    aux = c("sigma", "mix", "shiftprop"),
    aux_lb = c(0, 0, 0),
    aux_ub = c(Inf, 1, 1)
  ),
  "simplex" = list(
    constructor = simplex,
    rng = rsimplex,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf
  ),
  "softplusnormal" = list(
    constructor = softplusnormal,
    rng = rsoftplusnormal,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf,
    link = "softplus"
  ),
  "symlognormal" = list(
    constructor = symlognormal,
    rng = rsymlognormal,
    aux = "sigma",
    aux_lb = 0,
    aux_ub = Inf,
    link = "symlog"
  ),
  "unit_lindley" = list(
    constructor = unit_lindley,
    rng = runit_lindley,
    aux = character(0),
    aux_lb = numeric(0),
    aux_ub = numeric(0)
  ),
  "weibull" = list(
    constructor = function(link = NULL) brms::brmsfamily("weibull", link = link),
    rng = rweibull_median,
    aux = "shape",
    aux_lb = 0,
    aux_ub = Inf
  )
)

# Table of link and response functions used by link_lookup().
family_registry_links <- list(
  "logit" = list(link = logit, response = inv_logit),
  "cauchit" = list(link = cauchit, response = inv_cauchit),
  "cloglog" = list(link = cloglog, response = inv_cloglog),
  "identity" = list(link = identity, response = identity),
  "log" = list(link = log, response = exp),
  "softplus" = list(link = softplus, response = inv_softplus),
  "symlog" = list(link = symlog, response = inv_symlog)
)

# Fetch a registry entry, erroring clearly on unknown identifiers.
family_registry_entry <- function(family) {
  if (!is.character(family) || length(family) != 1L || is.na(family)) {
    stop("'family' must be a single string identifying a family.")
  }
  entry <- family_registry[[family]]
  if (is.null(entry)) {
    stop(
      "Unknown family identifier '", family, "'. Supported families: ",
      paste(names(family_registry), collapse = ", "),
      call. = FALSE
    )
  }
  entry
}

#' Lookup function for brms families via string identifier
#'
#' Looks up the family in the bayesfam family registry, which covers all
#' families exported by bayesfam plus a set of common brms built-in
#' families. Unknown identifiers raise an error.
#'
#' @include misc.R shifted_inv_gaussian.R shifted_lognormal_uniform.R simplex.R softplusnormal.R symlognormal.R unitlindley.R weibull_median.R
#'
#' @param family String identifier of the family.
#' @param link Link to be passed to the family function. If NULL (default),
#'   the constructor defaults are preserved.
#'
#' @return A brmsfamily object matching the string identifier and using the link
#' @export
#'
#' @examples
#' brms_family_lookup("weibull", "softplus")
brms_family_lookup <- function(family, link = NULL) {
  entry <- family_registry_entry(family)
  if (is.null(link)) {
    return(entry$constructor())
  }
  entry$constructor(link = link)
}

#' Lookup function for RNGs via string identifier
#'
#' Unknown identifiers raise an error.
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
  family_registry_entry(family)$rng
}

#' Lookup function for link and repsonse functions via string identifier.
#'
#' If a transformed normal likelihood is passed, the respective built-in link will
#' be returned instead of the identity link that would commonly be used with
#' transformed normal likelihood families.
#'
#' Unknown link or family identifiers raise an error.
#'
#' @param link String identifier for the link function of interest.
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
  if (!is.null(family)) {
    entry <- family_registry_entry(family)
    if (!is.null(entry$link)) {
      link <- entry$link
    }
  }
  link_functions <- family_registry_links[[link]]
  if (is.null(link_functions)) {
    stop(
      "Unknown link identifier '", link, "'. Supported links: ",
      paste(names(family_registry_links), collapse = ", "),
      call. = FALSE
    )
  }
  if (inv) {
    return(link_functions$response)
  }
  link_functions$link
}

#' Lookup function for the names of the auxiliary parameters of a likelihood
#'
#' Returns the auxiliary distributional parameter names of a family, i.e. all
#' dpars except mu. Zero-auxiliary families return an empty character vector.
#' Unknown identifiers raise an error.
#'
#' @param family The identifier string of a family.
#'
#' @return Character vector of auxiliary parameter names.
#' @export
#'
#' @examples
#' aux_family_parameters_lookup("beta")
aux_family_parameters_lookup <- function(family) {
  family_registry_entry(family)$aux
}

#' Lookup for limits of family auxiliary parameters.
#'
#' Returns a list with numeric lower and upper bound vectors of the auxiliary
#' parameters, matching the order of aux_family_parameters_lookup().
#' Zero-auxiliary families return empty numeric vectors. Unknown identifiers
#' raise an error.
#'
#' @param family The identifier string of a family.
#'
#' @return List containing lower and upper bounds for the auxiliary parameter.
#' @export
#'
#' @examples
#' aux_limits_lookup("beta")
aux_limits_lookup <- function(family) {
  entry <- family_registry_entry(family)
  list(lb = entry$aux_lb, ub = entry$aux_ub)
}
