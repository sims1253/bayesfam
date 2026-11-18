# Quantile function for the Inverse Burr distribution, with Median parametrization.

Quantile function for the Inverse Burr distribution, with Median
parametrization.

## Usage

``` r
qinverse_burr(p, mu, tau, gamma)
```

## Arguments

- p:

  Quantile to be calculated, p in (0, 1)

- mu:

  Median parameter, mu \> 0

- tau:

  Shape1 parameter, tau \> 0

- gamma:

  Shape2 parameter, gamma \> 0

## Value

Inverse of CDF, calculates a value, given a probability p

## Examples

``` r
p <- seq(from = 0.01, to = 0.99, length.out = 100)
plot(p, qinverse_burr(p, mu = 2, tau = 2, gamma = 3), type = "l")
```
