# brms family expect recovery. Tries linear baysian model y ~ 1. Checking of the arguments done in construct_brms.

brms family expect recovery. Tries linear baysian model y ~ 1. Checking
of the arguments done in construct_brms.

## Usage

``` r
expect_brms_family(
  n_data_sampels = 1000,
  intercept,
  ref_intercept = NULL,
  aux_par = NA,
  aux2_par = NA,
  rng_link,
  parameter_link,
  family,
  rng,
  aux_name = NULL,
  aux2_name = NULL,
  seed = 1235813,
  data_threshold = NULL,
  thresh = 0.05,
  debug = FALSE,
  formula = y ~ 1,
  prior = NULL
)
```

## Arguments

- n_data_sampels:

  How many samples per chain. Positive integer scalar. Default = 1000.

- intercept:

  Intercept for data generating RNG.

- ref_intercept:

  Reference intercept to compare model against. If NULL (default) uses
  the given intercept.

- aux_par:

  Auxiliary parameter of each distribution.

- aux2_par:

  2nd Auxiliary parameter (if applicable).

- rng_link:

  Link function pointer used for data generation. Mainly for transformed
  normal distributions.

- parameter_link:

  Link function pointer for the latent parameters. Used to transform for
  comparison with ref_intercept

- family:

  brms family under test.

- rng:

  function pointer of bespoke RNG for the family to be tested.

- aux_name:

  BRMS string of aux_par argument name. Single string.

- aux2_name:

  BRMS string of aux2_par argument name, if applicable. Single string.

- seed:

  Seed argument, so that input data is always the same in each test.
  brms test does not test RNG and is not guaranteed to fit on all data.
  Positive Integer scalar, Default = 1235813. Seed is stored before test
  and restored after it finished. If wants not to use a seed set to NA.

- data_threshold:

  Usually unused. But in rare cases, data too close at the boundary may
  cause trouble. If so, set a two entry real vector c(lower, upper). If
  one of them is NA, the data will not be capped for that boundary.
  Default = Null, will be in R terms "invisible" and will not cap any
  input data.

- thresh:

  Acceptable threshold for quantiles of recovered arguments. Scalar or
  2-entry real vector within (0, 1). Vector is used as is, scalar will
  be interpreted as c(thresh, 1-thresh). Default = 0.05

- debug:

  Scalar Boolean argument, whether debug info is printed or not. Default
  = False.

- formula:

  the formula used in the brms fit

- prior:

  any priors for brms

## Value

None

## Examples

``` r
result <- bayesfam:::expect_brms_family(
  intercept = 5,
  aux_par = 2,
  ref_intercept = 5,
  rng_link = identity,
  parameter_link = log,
  family = betaprime,
  rng = rbetaprime,
  aux_name = "phi"
)
#> Error in .fun(model_code = .x1) : 
#>   Boost not found; call install.packages('BH')
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
print(result)
#> Error: object 'result' not found
```
