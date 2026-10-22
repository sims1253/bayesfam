# Lognormal Natural density function with value space parameterisation (mu same as x)

Lognormal Natural density function with value space parameterisation (mu
same as x)

## Usage

``` r
dlognormal_natural(x, mu, sigma, log = FALSE)
```

## Arguments

- x:

  vector of values, x \> 0

- mu:

  Mean, mu \> 0

- sigma:

  Shape, sigma \> 0

- log:

  logical; if TRUE, log(pdf) is returned

## Value

Density of the Lognormal Natural given x, mu and sigma.

## Examples

``` r
x <- seq(from = 0.1, to = 5, length.out = 100)
plot(x, dlognormal_natural(x, 1, 2))
```
