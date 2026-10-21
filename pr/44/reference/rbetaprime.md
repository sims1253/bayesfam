# RNG for the beta prime distribution

RNG for the beta prime distribution

## Usage

``` r
rbetaprime(n, mu = 1, phi = 1)
```

## Arguments

- n:

  Number of samples.

- mu:

  Mean, mu \> 0.

- phi:

  Precision, phi \> 0.

## Value

Random numbers from the beta prime distribution.

## Examples

``` r
hist(rbetaprime(100, mu = 1, phi = 2))
```
