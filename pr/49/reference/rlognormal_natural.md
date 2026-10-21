# Lognormal Natural RNG function

Lognormal Natural RNG function

## Usage

``` r
rlognormal_natural(n, mu = 1, sigma = 1)
```

## Arguments

- n:

  number of observations

- mu:

  mean, mu \> 0

- sigma:

  sigma, sigma \> 0

## Value

n samples drawn from the Lognormal natural distribution

## Examples

``` r
hist(rlognormal_natural(100, 1, 2))
```
