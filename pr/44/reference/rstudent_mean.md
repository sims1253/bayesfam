# Custom rstudent with default df value mu and df arguments switched to comply to the Bayesfam allowing compatibility for Bayesim

Custom rstudent with default df value mu and df arguments switched to
comply to the Bayesfam allowing compatibility for Bayesim

## Usage

``` r
rstudent_mean(n, mu = 0, df = 1, sigma = 1)
```

## Arguments

- n:

  Number of sampels to draw, has to be a scalar natural

- mu:

  Mean argument, mu unbound

- df:

  Degrees of freedom variable

- sigma:

  Shape parameter, sigma \> 0

## Value

Vector of length n in student distribution

## Examples

``` r
hist(rstudent_mean(100, 1, 2, 1))
```
