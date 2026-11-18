# RNG for the gumbel distribution

RNG for the gumbel distribution

## Usage

``` r
rgumbel_mean(n, mu = 0, sigma = 1)
```

## Arguments

- n:

  Number of samples

- mu:

  Mean, unbound

- sigma:

  Scale, sigma \> 0

## Value

Random numbers from the gumbel distribution.

## Examples

``` r
hist(rgumbel_mean(100, mu = 2, sigma = 2))
```
