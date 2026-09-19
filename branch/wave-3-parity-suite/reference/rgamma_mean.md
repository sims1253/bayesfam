# Mean parameterization of the gamma RNG.

Mean parameterization of the gamma RNG.

## Usage

``` r
rgamma_mean(n, mu = 1, a = 1)
```

## Arguments

- n:

  Number of samples, scalar natural number.

- mu:

  Mean parameter, mu \> 0.

- a:

  Shape parameter, a \> 0.

## Value

n Gamma distributed samples.

## Examples

``` r
hist(log(rgamma_mean(10000, mu = 2, a = 1)))
```
