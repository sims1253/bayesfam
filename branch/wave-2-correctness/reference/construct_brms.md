# Construct brms family for simple linear y ~ 1 model.

Construct brms family for simple linear y ~ 1 model.

## Usage

``` r
construct_brms(
  n_data_sampels,
  intercept,
  aux_par = NA,
  aux2_par = NA,
  rng_link,
  family,
  rng,
  seed = NULL,
  data_threshold = NULL,
  formula = y ~ 1,
  prior = NULL
)
```

## Arguments

- n_data_sampels:

  How many samples per chain. Positive integer scalar.

- intercept:

  Intercept data argument, real scalar.

- aux_par:

  aux_par argument of each distribution.

- aux2_par:

  aux2_par argument applicable distribution

- rng_link:

  Link function pointer used data. For positive bounded uses exp as
  example.

- family:

  brms family under test.

- rng:

  function pointer of bespoke RNG for the family to be tested.

- seed:

  Seed argument, so that input data is always the same in each test.
  brms test does not test RNG and is not guaranteed to fit on all data.
  Positive Integer scalar, Default = NA will do nothing. Seed is stored
  before and restored after.

- data_threshold:

  Usually unused. But in rare cases, data too close at the boundary may
  cause trouble. If so, set a two entry real vector c(lower, upper). If
  one of them is NA, the data will not be capped for that boundary.
  Default = Null, will be in R terms "invisible" and will not cap any
  input data.

- formula:

  the formula used in the brms fit

- prior:

  any priors for brms

## Value

brms model for the specified family.

## Examples

``` r
posterior_fit <- bayesfam:::construct_brms(
  n_data_sampels = 1000,
  intercept = 5.0,
  aux_par = 2.0,
  rng_link = identity,
  family = betaprime,
  rng = rbetaprime
)
#> Error in .fun(model_code = .x1) : 
#>   Boost not found; call install.packages('BH')
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
plot(posterior_fit)
#> Error: object 'posterior_fit' not found
# beta_prime uses log-link for Intercept
```
