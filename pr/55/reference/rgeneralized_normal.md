# RNG for the generalized_normal distribution

RNG for the generalized_normal distribution

## Usage

``` r
rgeneralized_normal(n, mu = 0, sigma = 1, beta = 1)
```

## Arguments

- n:

  Number of sample

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

- beta:

  shape, beta \> 0

## Value

Random numbers from the generalized_normal distribution.

## Examples

``` r
hist(rgeneralized_normal(100, mu = 2, sigma = 2, beta = 2))
```
