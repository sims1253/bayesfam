# Lognormal RNG-function in median parametrization.

Lognormal RNG-function in median parametrization.

## Usage

``` r
rlognormal(n, mu = 0, sigma = 1)
```

## Arguments

- n:

  Number of draws

- mu:

  Median parameter, mu unbound, mu already log transformed

- sigma:

  Sigma shape parameter, sigma \> 0

## Value

n Lognormal distributed samples

## Examples

``` r
hist(rlognormal(100, 1, 0.5))
```
