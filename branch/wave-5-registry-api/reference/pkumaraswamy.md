# Kumaraswamy CDF in median parametrisation

Kumaraswamy CDF in median parametrisation

## Usage

``` r
pkumaraswamy(x, mu = 0.5, p = 1)
```

## Arguments

- x:

  CDF of x over lower tail, x e (0, 1)

- mu:

  Median parameter, mu e (0, 1)

- p:

  shape parameter, p \> 0

## Value

p(x \| mu, p)

## Examples

``` r
x <- seq(from = 0.01, to = 0.99, length.out = 1000)
plot(x, pkumaraswamy(x, mu = 0.5, p = 1), type = "l")
```
