# Test if an RNG asymptotically approaches the true location

Is currently not in a reliable state useful for testing. Needs more
thought.

## Usage

``` r
test_rng_asym(
  rng_fun,
  metric_mu,
  n_samples = c(10, 10000),
  mu_list,
  aux_list,
  aux2_list = NA,
  mu_link = identity,
  allowed_failures = 0.05
)
```

## Arguments

- rng_fun:

  RNG function under test

- metric_mu:

  Metric to be used on RNG data (usually mean or median)

- n_samples:

  Default=c(10, 10000), sample sizes for rng test, len \>= 2. Gets
  sorted in routine.

- mu_list:

  Metric data used as RNG argument and to be compared to (usually mean
  or median)

- aux_list:

  Auxiliary parameter value list.

- aux2_list:

  Auxiliary parameter for second parameter (if applicable).

- mu_link:

  Default=identity, optional link-function argument, for example

- allowed_failures:

  Default=0.05 marks, that 5 percent of all test cases are allowed to
  fail for the whole test to still succeed.

## Value

Success or failure with message

## Examples

``` r
eps <- 0.001
mu_list <- seq(from = 1 + eps, to = 20, length.out = 10)
phis <- seq(from = 2 + eps, to = 20, length.out = 10)
result <- bayesfam:::test_rng_asym(
  rng_fun = rbetaprime,
  metric_mu = mean,
  mu_list = mu_list,
  aux_list = phis,
)
print(result)
#> <expectation_success/expectation/condition>
#> As expected
```
