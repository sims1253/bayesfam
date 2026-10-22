# Contract fixtures for the custom brms families (#28).
#
# One entry per family, named like its `log_lik_<name>` callback: a response
# domain (support and three values inside it) and two valid values per
# distributional parameter, in the order of the family's dpars. Everything
# else (brms family, density, RNG, Stan code) is looked up by name.
# Optional fields: `family` overrides the constructor name, `vreal` holds
# per-observation data, `lpdf` and `rng` replace `d<family>(y, ..., log = TRUE)`
# and `r<family>(n, ...)`; all take the dpars and then the vreal data.
fixture <- function(domain, ...) c(domain, list(pars = list(...)))
unit <- list(support = c(0, 1), y = c(0.2, 0.5, 0.8))
positive <- list(support = c(0, Inf), y = c(0.5, 1, 3))
real <- list(support = c(-Inf, Inf), y = c(-1, 0.5, 2))

# fmt: skip
family_fixtures <- list(
  betaprime = fixture(positive, mu = c(1, 2), phi = c(2, 5)),
  cauchitnormal = fixture(unit, mu = c(-0.5, 0.5), sigma = c(0.5, 1)),
  cloglognormal = fixture(unit, mu = c(-1, -0.5), sigma = c(0.3, 0.5)),
  generalized_gamma = fixture(positive, mu = c(0, 0.5), sigma = c(0.5, 1), Q = c(-0.5, 1.5)),
  generalized_normal = fixture(real, mu = c(0, 1), sigma = c(1, 2), beta = c(1.5, 3)),
  gompertz = fixture(positive, mu = c(1, 2), beta = c(0.5, 1)),
  gumbel_mean = fixture(real, mu = c(0, 1), sigma = c(0.5, 1)),
  inverse_burr = fixture(positive, mu = c(1, 2), tau = c(2, 3), gamma = c(3, 4)),
  kumaraswamy = fixture(unit, mu = c(0.3, 0.6), p = c(2, 4)),
  logistic = fixture(real, mu = c(0, 1), sigma = c(0.5, 1)),
  logitnormal = fixture(unit, mu = c(-0.5, 0.5), sigma = c(0.5, 1)),
  lognormal_natural = fixture(positive, mu = c(1, 2), sigma = c(0.5, 1)),
  lomax = fixture(positive, mu = c(1, 2), alpha = c(3, 5)),
  shifted_inv_gaussian = fixture(list(support = c(0.2, Inf), y = c(0.5, 1, 3)), mu = c(1, 2), shape = c(1, 3), ndt = c(0.2, 0.2)),
  shifted_lognormal_uniform = c(
    fixture(list(support = c(0, Inf), y = c(0.05, 1, 12)), mu = c(0, 0.5), sigma = c(0.5, 1), mix = c(0.1, 0.3), shiftprop = c(0.2, 0.8)),
    # the d/r functions take the shift itself, and d always returns the log density
    vreal = list(list(vreal1 = 0.3, vreal2 = 10)),
    lpdf = function(y, mu, sigma, mix, shiftprop, max_shift, max_uniform) dshifted_lognormal_uniform(y, mu, sigma, mix, shiftprop * max_shift, max_uniform),
    rng = function(n, mu, sigma, mix, shiftprop, max_shift, max_uniform) rshifted_lognormal_uniform(n, mu, sigma, mix, shiftprop * max_shift, max_uniform)
  ),
  simplex = fixture(unit, mu = c(0.3, 0.6), sigma = c(0.5, 1)),
  softplusnormal = fixture(positive, mu = c(0.5, 1), sigma = c(0.5, 1)),
  symlognormal = fixture(real, mu = c(0, 1), sigma = c(0.5, 1)),
  unitlindley = c(fixture(unit, mu = c(0.3, 0.6)), family = "unit_lindley")
)

# The fixture of family `f`, completed with its brms family, lpdf and rng.
family_fixture <- function(f) {
  fx <- family_fixtures[[f]]
  ns <- asNamespace("bayesfam")
  name <- if (is.null(fx$family)) f else fx$family
  fx$fam <- get(name, ns)()
  dens <- get(paste0("d", name), ns)
  if (is.null(fx$lpdf)) {
    fx$lpdf <- function(y, ...) dens(y, ..., log = TRUE)
  }
  if (is.null(fx$rng)) {
    fx$rng <- get(paste0("r", name), ns)
  }
  fx
}
