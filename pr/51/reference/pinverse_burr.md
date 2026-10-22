# Cumulative distribution function for the Inverse Burr distribution, with Median parametrization.

Cumulative distribution function for the Inverse Burr distribution, with
Median parametrization.

## Usage

``` r
pinverse_burr(q, mu, tau, gamma)
```

## Arguments

- q:

  Value to evaluate the CDF at, q \> 0

- mu:

  Median parameter, mu \> 0

- tau:

  Shape1 parameter, tau \> 0

- gamma:

  Shape2 parameter, gamma \> 0

## Value

CDF of Inverse Burr Distribution

## Examples

``` r
x <- seq(from = 0.1, to = 10, length.out = 100)
plot(x, pinverse_burr(x, mu = 2, tau = 2, gamma = 3), type = "l")
```
