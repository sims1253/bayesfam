# Simplex RNG function in Median parametrization.

Based on code from simplexreg Peng Zhang, Zhenguo Qiu, Chengchun Shi
(2016). simplexreg: An R Package for Regression Analysis of Proportional
Data Using the Simplex Distribution. Journal of Statistical Software,
71(11), 1-21. doi:10.18637/jss.v071.i11

## Usage

``` r
rsimplex(n, mu = 0.5, sigma = 1)
```

## Arguments

- n:

  Number of samples to draw, as a natural number scalar.

- mu:

  Mean parameter, mu e (0, 1)

- sigma:

  shape parameter, Sigma unbound

## Value

n samples in Simplex distribution.

## Examples

``` r
hist(rsimplex(10000, mu = 0.7, sigma = 2))
```
