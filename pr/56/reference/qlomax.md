# Quantile function for the Lomax distribution, with Mean parametrization.

Quantile function for the Lomax distribution, with Mean parametrization.

## Usage

``` r
qlomax(p, mu, alpha)
```

## Arguments

- p:

  Quantile to be calculated

- mu:

  Median argument of Lomax

- alpha:

  Alpha argument of Gompertz

## Value

Inverse of CDF, calculates a value, given a probability p

## Examples

``` r
x <- seq(from = 0, to = 1, length.out = 100)
plot(x, qlomax(x, mu = 1, alpha = 2), type = "l")
```
