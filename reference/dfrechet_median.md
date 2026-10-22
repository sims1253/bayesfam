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

  Shape, nu \> 0

## Value

dfrechet(x \| mu, nu)

## Details

Define scale parameter sigma as \$\$\sigma(\mu, \nu) := \mu \cdot
\log(2)^{1 / \nu}\$\$

The Frechet distribution has density \$\$f(y) = (\nu /\sigma) \* (y /
\sigma)^{-(\nu + 1)} \* exp(-(y / \sigma)^{-\nu}) \$\$

## Examples

``` r
x <- seq(from = 0.1, to = 20, length.out = 1000)
plot(x, dfrechet_median(x, mu = 6, nu = 4), type = "l")
```
