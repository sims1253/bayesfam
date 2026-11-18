# RNG function for the Gompertz distribution, with Median parametrization.

RNG function for the Gompertz distribution, with Median parametrization.

## Usage

``` r
rgompertz(n, mu = 1, beta = 0.5)
```

## Arguments

- n:

  Number of draws

- mu:

  Median parameter

- beta:

  Scale parameter

## Value

A Gompertz distributed RNG vector of size n

## Examples

``` r
hist(rgompertz(n = 100, mu = 2, beta = 0.1))
```
