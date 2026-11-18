# Median parameterization of the Weibull RNG.

Median parameterization of the Weibull RNG.

## Usage

``` r
rweibull_median(n, mu = 1, k = 1)
```

## Arguments

- n:

  Number of samples, scalar natural number.

- mu:

  Median parameter, mu \> 0.

- k:

  Shape parameter, k \> 0.

## Value

n Weibull distributed samples.

## Examples

``` r
hist(log(rweibull_median(10000, mu = 2, k = 1)))
```
