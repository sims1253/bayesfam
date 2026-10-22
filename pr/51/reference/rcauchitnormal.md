# Cauchitnormal RNG-function

Cauchitnormal RNG-function

## Usage

``` r
rcauchitnormal(n, mu = 0, sigma = 1)
```

## Arguments

- n:

  Number of draws

- mu:

  Median parameter, mu unbound, mu already cauchit transformed

- sigma:

  Sigma shape parameter, sigma \> 0

## Value

n chauchitnormal distributed samples

## Examples

``` r
hist(rcauchitnormal(100, 0.5, 2))
```
