# Shifted inverse Gaussian RNG function

Shifted inverse Gaussian RNG function

## Usage

``` r
rshifted_inv_gaussian(n, mu = 1, shape = 1, shift = 1)
```

## Arguments

- n:

  number of observations

- mu:

  Mean, mu \> 0

- shape:

  Shape, shape unbound

- shift:

  Shift, shift \>= 0

## Value

n samples drawn from the shifted inverse Gaussian distribution

## Examples

``` r
hist(rshifted_inv_gaussian(100, 1, 1, 1))
```
