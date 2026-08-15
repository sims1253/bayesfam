# Custom rexgauss with default values

Custom rexgauss with default values

## Usage

``` r
rexgauss_mean(n, mu = 0, sigma = 1, beta = 1)
```

## Arguments

- n:

  Number of sampels to draw, has to be a scalar natural

- mu:

  Mean argument, mu unbound

- sigma:

  Shape parameter, sigma \> 0

- beta:

  Shape parameter, beta \> 0

## Value

Vector of length n in exgaussian distribution

## Examples

``` r
hist(rexgauss_mean(100, 1, 2, 1))
```
