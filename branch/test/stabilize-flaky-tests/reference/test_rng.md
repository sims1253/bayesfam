# Tests if an RNG can achieve a good enough location often enough

Tests if an RNG can achieve a good enough location often enough

## Usage

``` r
test_rng(
  rng_fun,
  metric_mu,
  n,
  mu_list,
  aux_list = NA,
  aux2_list = NA,
  mu_eps,
  p_acceptable_failures,
  mu_link = identity,
  relative = FALSE,
  debug = TRUE
)
```

## Arguments

- rng_fun:

  RNG function under test

- metric_mu:

  Metric to be used on RNG data (usually mean or median)

- n:

  Sample size for the rng test.

- mu_list:

  Metric data used as RNG argument and to be compared to (usually mean
  or median)

- aux_list:

  Auxiliary parameter

- aux2_list:

  Auxiliary parameter for second parameter (if applicable).

- mu_eps:

  Acceptable difference of \|mu - metric_mu(rng_fun)

- p_acceptable_failures:

  Acceptable rate of failure, relative value of difference bigger mu_eps

- mu_link:

  Default=identity, optional link-function argument, for example useful
  in link-normal-distributions

- relative:

  True if the error should be relative to the mu_list, Default = FALSE

- debug:

  bool argument, if set will printout the difference, in case of
  failure. Default = FALSE

## Value

Nothing actually, just wraps the test

## Examples

``` r
eps <- 1e-6
mu_list <- seq(from = 1 + eps, to = 20, length.out = 10)
phis <- seq(from = 2 + eps, to = 20, length.out = 10)
result <- bayesfam:::test_rng(
  rng_fun = rbetaprime, metric_mu = mean, n = 10000, mu_list = mu_list,
  aux_list = phis, mu_eps = 0.2, p_acceptable_failures = 0.05
)
print(result)
#> <expectation_success/expectation/condition>
#> As expected
```
