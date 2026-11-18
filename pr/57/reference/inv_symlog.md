# Symlog response function

Symlog response function

## Usage

``` r
inv_symlog(x)
```

## Source

Based on Hafner, D., Pasukonis, J., Ba, J., & Lillicrap, T. (2023).
Mastering Diverse Domains through World Models.
(<https://doi.org/10.48550/arXiv.2301.04104>)

## Arguments

- x:

  value to be transformed, x is unbound

## Value

inv_symlog of x, result is unbound

## Examples

``` r
inv_symlog(0)
#> [1] 0
inv_symlog(10)
#> [1] 22025.47
inv_symlog(-10)
#> [1] -22025.47
```
