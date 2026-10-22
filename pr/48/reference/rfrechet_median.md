# Median parameterization of the Fréchet RNG

Median parameterization of the Fréchet RNG

## Usage

``` r
rfrechet_median(n, mu = 1, nu = 2)
```

## Arguments

- n:

  Number samples to draw

- mu:

  Median

- nu:

  Shape, nu \> 0

## Value

n samples in Frechet-Distribution

## Examples

``` r
hist(rfrechet_median(100, mu = 1, nu = 2))
```
