# RNG function for the Inverse Burr distribution, with Median parametrization.

RNG function for the Inverse Burr distribution, with Median
parametrization.

## Usage

``` r
rinverse_burr(n, mu = 1, tau = 1, gamma = 1)
```

## Source

Parameterization follows the actuar package
(<https://search.r-project.org/CRAN/refmans/actuar/html/InverseBurr.html>).

## Arguments

- n:

  Number of draws

- mu:

  Median parameter, mu \> 0

- tau:

  Shape1 parameter, tau \> 0 (actuar: shape1)

- gamma:

  Shape2 parameter, gamma \> 0 (actuar: shape2)

## Value

An Inverse Burr distributed RNG vector of size n

## Examples

``` r
hist(rinverse_burr(1000, mu = 2, tau = 2, gamma = 3))
```
