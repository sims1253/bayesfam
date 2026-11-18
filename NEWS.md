# bayesfam (development version)

## Breaking changes

* `lomax()` now uses `link_alpha = "logm1"` by default; the old default was
  `"log1p"`. The old link let `alpha` fall below 1, where the Lomax density
  is undefined. Code that calls `lomax()` without `link_alpha` now fits
  `alpha` on a different link scale, which changes its priors and
  coefficients; pass `link_alpha = "log1p"` to get the old behaviour.

## Bug fixes

* `rlomax()` now rejects `alpha <= 1`, like `dlomax()` and `qlomax()`.
* `cauchitnormal()`, `cloglognormal()`, `logitnormal()` and
  `softplusnormal()` no longer bound `mu`. It is the location of the normal
  variable before the transform, so it can take any real value.
* `dgompertz()`, `qgompertz()`, `dsoftplusnormal()` and `rsoftplusnormal()`
  now reject a vector argument when any element of it is invalid, not only
  when all elements are.
* `dgompertz()`, `qgompertz()`, `rgompertz()` and the Gompertz Stan code use
  `expm1()`, so they stay accurate when `mu * beta` is small; before,
  `qgompertz(0.5, 1, 1e-16)` gave 0 instead of the median 1.
* `runit_lindley()` and the unit Lindley Stan RNG draw from a gamma
  mixture; for `mu` below about 0.0013 every draw was `NaN`.
* `rshifted_inv_gaussian()`, `rsimplex()` and their Stan RNGs compute the
  smaller inverse Gaussian root without cancellation; for small shapes or
  large `sigma`, draws landed at or below the shift, or at 0.
* `inv_cloglog()` uses `expm1()`, so it no longer rounds to 0 below
  `x` of about -37.
* `lomax()` includes brms' `logm1` Stan functions only when `link` or
  `link_alpha` is `"logm1"`.
* `cloglog()` and `cauchit()` give their own name in their error messages;
  before, both said "logit link".
* The unit Lindley `posterior_predict` reads `mu` through `brms::get_dpar()`.

## New features

* `posterior_epred` is now implemented for `gompertz()`,
  `generalized_gamma()`, `cauchitnormal()`, `cloglognormal()` and
  `softplusnormal()`. The generalized gamma mean is `Inf` where it does not
  exist, for `Q * sigma <= -1`, and stays accurate as `Q` approaches 0, down
  to the smallest normal double; for subnormal `Q` it overflows to `Inf`. The
  Gompertz mean stays finite where `expm1(mu * beta)` overflows, as long as
  the product `mu * beta` is itself finite.
* `unit_lindley()`, `generalized_gamma()`, `simplex()`,
  `shifted_inv_gaussian()` and `shifted_lognormal_uniform()` now define a
  Stan `_rng` function.
* `dfrechet_median()` gains a `log` argument.
