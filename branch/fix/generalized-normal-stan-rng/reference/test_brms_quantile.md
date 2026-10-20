# Check, that data of the posterior is close enough to the reference data.

Check, that data of the posterior is close enough to the reference data.

## Usage

``` r
test_brms_quantile(posterior_data, arg_name, reference, thresh, debug = FALSE)
```

## Arguments

- posterior_data:

  Data fitted and drawn, a brms object, which is a list in R terms

- arg_name:

  Name of the argument variable to check, as single string

- reference:

  Reference value to check against, single real scalar

- thresh:

  real scalar or 2-length vector of quantile bounds. For scalar
  constructs bound as `[thresh, 1-thresh]` thresh has to be inside the
  Unit-Interval.

- debug:

  True for verbose output of test results.

## Value

Single boolean success, fail or error

## Examples

``` r
fit <- bayesfam:::construct_brms(
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
result <- bayesfam:::test_brms_quantile(
  posterior_data = fit, arg_name = "phi", 2.0, 0.025
)
#> Error: object 'fit' not found
plot(fit)
#> Error: object 'fit' not found
# beta_prime uses log-link for Intercept
```
