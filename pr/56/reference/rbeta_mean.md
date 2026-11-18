# Beta distribution RNG for mean parameterization

Beta distribution RNG for mean parameterization

## Usage

``` r
rbeta_mean(n, mu = 0.5, phi = 4)
```

## Arguments

- n:

  Number of draws.

- mu:

  Mean

- phi:

  Precision

## Value

n samples Beta distributed.

## Examples

``` r
hist(rbeta_mean(1000, mu = 0.5, phi = 1))
```
