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
* `dgompertz()` uses `expm1()`, so it stays accurate when `mu * beta` is
  small.
* `cloglog()` and `cauchit()` give their own name in their error messages;
  before, both said "logit link".
* The unit Lindley `posterior_predict` reads `mu` through `brms::get_dpar()`.

## New features

* `posterior_epred` is now implemented for `gompertz()`,
  `generalized_gamma()`, `cauchitnormal()`, `cloglognormal()` and
  `softplusnormal()`. The generalized gamma mean is `Inf` where it does not
  exist, for `Q <= -1 / sigma`.
* `unit_lindley()`, `generalized_gamma()`, `simplex()`,
  `shifted_inv_gaussian()` and `shifted_lognormal_uniform()` now define a
  Stan `_rng` function.
* `dfrechet_median()` gains a `log` argument.
