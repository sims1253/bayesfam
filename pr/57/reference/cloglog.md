# Complementary-Log-Log link function

Complementary-Log-Log link function

## Usage

``` r
cloglog(x)
```

## Arguments

- x:

  value of x to be transformed, x e (0, 1)

## Value

cloglog value of x, result is unbound

## Examples

``` r
x <- seq(from = 0.1, to = 0.9, length.out = 100)
plot(x, cloglog(x), type = "l")
```
