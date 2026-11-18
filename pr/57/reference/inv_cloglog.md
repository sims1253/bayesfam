# Complementary-Log-Log response function

Complementary-Log-Log response function

## Usage

``` r
inv_cloglog(x)
```

## Arguments

- x:

  value of x to be transformed, x unbound

## Value

inverse-cloglog value of x, result is e (0, 1)

## Examples

``` r
x <- seq(from = -3, to = 1, length.out = 100)
plot(x, inv_cloglog(x), type = "l")
```
