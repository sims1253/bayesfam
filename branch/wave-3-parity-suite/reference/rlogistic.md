# RNG for the logistic distribution

RNG for the logistic distribution

## Usage

``` r
rlogistic(n, mu = 0, sigma = 1)
```

## Arguments

- n:

  Number of samples.

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

## Value

Random numbers from the logistic distribution.

## Examples

``` r
hist(rlogistic(100, mu = 2, sigma = 2))
```
