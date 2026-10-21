# Median parameterization of the Fréchet pdf.

Median parameterization of the Fréchet pdf.

## Usage

``` r
dfrechet_median(x, mu, nu)
```

## Arguments

- x:

  x value space, x \> 0

- mu:

  Median

- nu:

  Shape

## Value

dfrechet(x \| mu, nu)

## Details

Define scale parameter sigma as \$\$\sigma(\mu, \nu) := \mu / \Gamma(1 -
1 / \nu)\$\$

The Frechet distribution has density \$\$f(y) = (\nu /\sigma) \* (y /
\sigma)^{-(1 - \nu)} \* exp(-(y / \sigma)^{-\nu}) \$\$

## Examples

``` r
x <- seq(from = 0.1, to = 20, length.out = 1000)
plot(x, dfrechet_median(x, mu = 6, nu = 4), type = "l")
```
