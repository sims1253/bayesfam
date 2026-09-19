# Symlog link function

Symlog link function

## Usage

``` r
symlog(x)
```

## Source

Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
Mastering Diverse Domains through World Models.
(<https://doi.org/10.48550/arXiv.2301.04104>)

## Arguments

- x:

  value to be transformed, x is unbound

## Value

symlog of x, result is unbound

## Examples

``` r
symlog(0)
#> [1] 0
symlog(1e10)
#> [1] 23.02585
symlog(-1e10)
#> [1] -23.02585
```
