# Softplus RNG-function in median parametrization.

Softplus RNG-function in median parametrization.

## Usage

``` r
rsoftplusnormal(n, mu = 1, sigma = 1)
```

## Arguments

- n:

  Number of draws

- mu:

  Median parameter, mu unbound, mu already log transformed

- sigma:

  Sigma shape parameter, sigma \> 0

## Value

n Softplus distributed samples

## Examples

``` r
hist(rsoftplusnormal(100, 1, 2))
```
