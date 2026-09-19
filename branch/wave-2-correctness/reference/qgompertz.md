# Quantile function for the Gompertz distribution, with Median parametrization.

Quantile function for the Gompertz distribution, with Median
parametrization.

## Usage

``` r
qgompertz(p, mu, beta)
```

## Arguments

- p:

  Quantile to be calculated

- mu:

  Median parameter

- beta:

  Scale parameter

## Value

Inverse of CDF, calculates a value, given a probability p

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 100)
plot(x, qgompertz(x, mu = 10, beta = 1), type = "l")
```
