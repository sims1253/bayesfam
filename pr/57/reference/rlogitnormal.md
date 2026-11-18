# Logitnormal RNG-function in median parametrization.

Logitnormal RNG-function in median parametrization.

## Usage

``` r
rlogitnormal(n, mu = 0, sigma = 1)
```

## Arguments

- n:

  number of observations

- mu:

  Median parameter, mu unbound, mu already logit transformed

- sigma:

  Sigma shape parameter, sigma \> 0

## Value

n Logitnormal distributed samples

## Examples

``` r
hist(rlogitnormal(100, 0.5, 2))
```
