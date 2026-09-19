# RNG function for the Lomax distribution, with Mean parametrization.

RNG function for the Lomax distribution, with Mean parametrization.

## Usage

``` r
rlomax(n, mu = 1, alpha = 10)
```

## Arguments

- n:

  Number of draws

- mu:

  Median argument of Lomax

- alpha:

  Eta argument of Lomax

## Value

A Lomax distributed RNG vector of size n

## Examples

``` r
hist(log(rlomax(1000, mu = 1, alpha = 2)))
```
