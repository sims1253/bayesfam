# Changelog

## bayesfam (development version)

### Breaking changes

- [`lomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/lomax.md)
  now uses `link_alpha = "logm1"` by default; the old default was
  `"log1p"`. The old link let `alpha` fall below 1, where the Lomax
  density is undefined. Code that calls
  [`lomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/lomax.md)
  without `link_alpha` now fits `alpha` on a different link scale, which
  changes its priors and coefficients; pass `link_alpha = "log1p"` to
  get the old behaviour.

### Bug fixes

- [`rlomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/rlomax.md)
  now rejects `alpha <= 1`, like
  [`dlomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/dlomax.md)
  and
  [`qlomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/qlomax.md).
- [`cauchitnormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/cauchitnormal.md),
  [`cloglognormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/cloglognormal.md),
  [`logitnormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/logitnormal.md)
  and
  [`softplusnormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/softplusnormal.md)
  no longer bound `mu`. It is the location of the normal variable before
  the transform, so it can take any real value.
- [`dgompertz()`](https://sims1253.github.io/bayesfam/pr/56/reference/dgompertz.md),
  [`qgompertz()`](https://sims1253.github.io/bayesfam/pr/56/reference/qgompertz.md),
  [`dsoftplusnormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/dsoftplusnormal.md)
  and
  [`rsoftplusnormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/rsoftplusnormal.md)
  now reject a vector argument when any element of it is invalid, not
  only when all elements are; so do
  [`dlomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/dlomax.md)
  and
  [`qlomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/qlomax.md)
  for `alpha`.
- [`dgompertz()`](https://sims1253.github.io/bayesfam/pr/56/reference/dgompertz.md),
  [`qgompertz()`](https://sims1253.github.io/bayesfam/pr/56/reference/qgompertz.md),
  [`rgompertz()`](https://sims1253.github.io/bayesfam/pr/56/reference/rgompertz.md)
  and the Gompertz Stan code use
  [`expm1()`](https://rdrr.io/r/base/Log.html), so they stay accurate
  when `mu * beta` is small; before, `qgompertz(0.5, 1, 1e-16)` gave 0
  instead of the median 1.
- [`runit_lindley()`](https://sims1253.github.io/bayesfam/pr/56/reference/runit_lindley.md)
  and the unit Lindley Stan RNG draw from a gamma mixture; for `mu`
  below about 0.0013 every draw was `NaN`.
- [`rshifted_inv_gaussian()`](https://sims1253.github.io/bayesfam/pr/56/reference/rshifted_inv_gaussian.md),
  [`rsimplex()`](https://sims1253.github.io/bayesfam/pr/56/reference/rsimplex.md)
  and their Stan RNGs compute the smaller inverse Gaussian root without
  cancellation; for small shapes or large `sigma`, draws landed at or
  below the shift, or at 0.
- [`inv_cloglog()`](https://sims1253.github.io/bayesfam/pr/56/reference/inv_cloglog.md)
  uses [`expm1()`](https://rdrr.io/r/base/Log.html), so it no longer
  rounds to 0 below `x` of about -37.
- [`lomax()`](https://sims1253.github.io/bayesfam/pr/56/reference/lomax.md)
  includes brms’ `logm1` Stan functions only when `link` or `link_alpha`
  is `"logm1"`.
- [`cloglog()`](https://sims1253.github.io/bayesfam/pr/56/reference/cloglog.md)
  and
  [`cauchit()`](https://sims1253.github.io/bayesfam/pr/56/reference/cauchit.md)
  give their own name in their error messages; before, both said “logit
  link”.
- The unit Lindley `posterior_predict` reads `mu` through
  [`brms::get_dpar()`](https://paulbuerkner.com/brms/reference/get_dpar.html).

### New features

- `posterior_epred` is now implemented for
  [`gompertz()`](https://sims1253.github.io/bayesfam/pr/56/reference/gompertz.md),
  [`generalized_gamma()`](https://sims1253.github.io/bayesfam/pr/56/reference/generalized_gamma.md),
  [`cauchitnormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/cauchitnormal.md),
  [`cloglognormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/cloglognormal.md)
  and
  [`softplusnormal()`](https://sims1253.github.io/bayesfam/pr/56/reference/softplusnormal.md).
  The generalized gamma mean is `Inf` where it does not exist, for
  `Q * sigma <= -1`, and stays accurate as `Q` approaches 0, down to the
  smallest normal double; for subnormal `Q` it overflows to `Inf`. The
  Gompertz mean stays finite where `expm1(mu * beta)` overflows, as long
  as the product `mu * beta` is itself finite.
- [`unit_lindley()`](https://sims1253.github.io/bayesfam/pr/56/reference/unit_lindley.md),
  [`generalized_gamma()`](https://sims1253.github.io/bayesfam/pr/56/reference/generalized_gamma.md),
  [`simplex()`](https://sims1253.github.io/bayesfam/pr/56/reference/simplex.md),
  [`shifted_inv_gaussian()`](https://sims1253.github.io/bayesfam/pr/56/reference/shifted_inv_gaussian.md)
  and
  [`shifted_lognormal_uniform()`](https://sims1253.github.io/bayesfam/pr/56/reference/shifted_lognormal_uniform.md)
  now define a Stan `_rng` function.
- [`dfrechet_median()`](https://sims1253.github.io/bayesfam/pr/56/reference/dfrechet_median.md)
  gains a `log` argument.
