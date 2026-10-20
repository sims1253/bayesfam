# Cloglognormal RNG-function

Cloglognormal RNG-function

## Usage

``` r
rcloglognormal(n, mu = -0.36, sigma = 0.75)
```

## Arguments

- n:

  Number of draws

- mu:

  Median parameter, mu unbound, mu already cloglog transformed

- sigma:

  Shape parameter

## Value

n cloglog-normally distributed samples

## Examples

``` r
hist(rcloglognormal(100, 0.5, 2))
```
