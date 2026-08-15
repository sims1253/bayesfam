# Kumaraswamy RNG function in Median parametrization.

Kumaraswamy RNG function in Median parametrization.

## Usage

``` r
rkumaraswamy(n, mu = 0.5, p = 2)
```

## Arguments

- n:

  number of observations

- mu:

  Median parameter, mu e (0, 1)

- p:

  Phi shape parameter, Phi \> 0

## Value

n samples in Kumaraswamy distribution.

## Examples

``` r
hist(rkumaraswamy(10000, mu = 0.5, p = 4))
```
